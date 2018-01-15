// Down-and-out European Continuous Barrier Option Pricing Code with reference to "down_and_out_adjusted_dynamic_prog.cpp"
// Written by Shihua Huang
#include <iostream>
#include <iomanip>
#include <cmath>
#include <fstream>
#include <cstdlib>
#include <algorithm>
#include <chrono>
#include <random>
#include "normdist.h"
#include "newmat.h"

using namespace std;

double up_factor, uptick_prob, risk_free_rate, strike_price;
double initial_stock_price, expiration_time, volatility, R, barrier_price;
int no_of_trials, no_of_discrete_barriers;
double simulate_call_option_price = 0.0;
double simulate_put_option_price = 0.0;
double simulate_call_option_price_adjusted = 0.0;
double simulate_put_option_price_adjusted = 0.0;

double max(double a, double b) {
	return (b < a) ? a : b;
}

std::default_random_engine generator;
double get_uniform()
{
	std::uniform_real_distribution <double> distribution(0.0, 1.0);
	double number = distribution(generator);
	return (number);
}

double N(const double& z) {
	if (z > 6.0) { return 1.0; }; // this guards against overflow 
	if (z < -6.0) { return 0.0; };
	double b1 = 0.31938153;
	double b2 = -0.356563782;
	double b3 = 1.781477937;
	double b4 = -1.821255978;
	double b5 = 1.330274429;
	double p = 0.2316419;
	double c2 = 0.3989423;
	double a = fabs(z);
	double t = 1.0 / (1.0 + a*p);
	double b = c2*exp((-z)*(z / 2.0));
	double n = ((((b5*t + b4)*t + b3)*t + b2)*t + b1)*t;
	n = 1.0 - b*n;
	if (z < 0.0) n = 1.0 - n;
	return n;
};

void adjusted_simulation()
{
	double delta_T = expiration_time / ((double)no_of_discrete_barriers);
	double delta_R = (risk_free_rate - 0.5*pow(volatility, 2))*delta_T;
	double delta_SD = volatility*sqrt(delta_T);
	//create dynamic arraies to store the value of mean and 
	//variance at each barrier point and calculate the probability for discrete pricing
	double *mean_at_sampling_instant1 = new double[no_of_discrete_barriers];
	double *mean_at_sampling_instant2 = new double[no_of_discrete_barriers];
	double *mean_at_sampling_instant3 = new double[no_of_discrete_barriers];
	double *mean_at_sampling_instant4 = new double[no_of_discrete_barriers];
	double *variance_at_sampling_instant = new double[no_of_discrete_barriers];

	for (int j = 0; j < no_of_trials; j++)
	{
		//create 4 variables for recording whether prices have breached barriers. 
		//They are breach_value1,2,3,4 and breach_value_adjusted1,2,3,4 for explicit 
		//simulation and (1-p)adjustment price respectively
		double current_stock_price1 = initial_stock_price;
		double current_stock_price2 = initial_stock_price;
		double current_stock_price3 = initial_stock_price;
		double current_stock_price4 = initial_stock_price;
		double breach_value1 = 1.0;
		double breach_value2 = 1.0;
		double breach_value3 = 1.0;
		double breach_value4 = 1.0;
		double breach_value_adjusted1 = 1.0;
		double breach_value_adjusted2 = 1.0;
		double breach_value_adjusted3 = 1.0;
		double breach_value_adjusted4 = 1.0;
		for (int i = 0; i < no_of_discrete_barriers; i++)
		{
			//use box-muller methods to simulate 4 price paths
			double x = get_uniform();
			double y = get_uniform();
			double a = sqrt(-2.0*log(x)) * cos(6.283185307999998*y);
			double b = sqrt(-2.0*log(x)) * sin(6.283185307999998*y);

			current_stock_price1 = current_stock_price1*exp(delta_R + delta_SD*a);
			current_stock_price2 = current_stock_price2*exp(delta_R - delta_SD*a);
			current_stock_price3 = current_stock_price3*exp(delta_R + delta_SD*b);
			current_stock_price4 = current_stock_price4*exp(delta_R - delta_SD*b);
			//see if stock price are not greter than barrier price, 
			//it is knocked out and record the breach_value as 0
			if (current_stock_price1 <= barrier_price)
				breach_value1 = 0.0;
			if (current_stock_price2 <= barrier_price)
				breach_value2 = 0.0;
			if (current_stock_price3 <= barrier_price)
				breach_value3 = 0.0;
			if (current_stock_price4 <= barrier_price)
				breach_value4 = 0.0;
		}
		//calculate the price in each trial as multiplying the payoff by 0 or 1 and divide 4
			simulate_call_option_price += (breach_value1*max(0.0, current_stock_price1 - strike_price) +
				breach_value2*max(0.0, current_stock_price2 - strike_price) +
				breach_value3*max(0.0, current_stock_price3 - strike_price) +
				breach_value4*max(0.0, current_stock_price4 - strike_price)) / 4.0;
			simulate_put_option_price += (breach_value1*max(0.0, strike_price - current_stock_price1) +
				breach_value2*max(0.0, strike_price - current_stock_price2) +
				breach_value3*max(0.0, strike_price - current_stock_price3) +
				breach_value4*max(0.0, strike_price - current_stock_price4)) / 4.0;
			//after simulating the 4-prices-path, we test the last price with barrier price,
			//if it breach the barrier price, the variable breach_value_adjusted should be recorded as 0
			if (current_stock_price1 <= barrier_price)
				breach_value_adjusted1 = 0.0;
			if (current_stock_price2 <= barrier_price)
				breach_value_adjusted2 = 0.0;
			if (current_stock_price3 <= barrier_price)
				breach_value_adjusted3 = 0.0;
			if (current_stock_price4 <= barrier_price)
				breach_value_adjusted4 = 0.0;
			//calculate the probability of breaching barrier price according to equation(2) in assignment instruction
			double probability_stock_path_has_not_hit_barrier_at_some_discrete_sampling_instant1 = 1.0;
			double probability_stock_path_has_not_hit_barrier_at_some_discrete_sampling_instant2 = 1.0;
			double probability_stock_path_has_not_hit_barrier_at_some_discrete_sampling_instant3 = 1.0;
			double probability_stock_path_has_not_hit_barrier_at_some_discrete_sampling_instant4 = 1.0;
			for (int i = 1; i <= no_of_discrete_barriers; i++)
			{
				mean_at_sampling_instant1[i] = initial_stock_price +
					(((double)i) / ((double)no_of_discrete_barriers)*(current_stock_price1 - initial_stock_price));
				mean_at_sampling_instant2[i] = initial_stock_price +
					(((double)i) / ((double)no_of_discrete_barriers)*(current_stock_price2 - initial_stock_price));
				mean_at_sampling_instant3[i] = initial_stock_price +
					(((double)i) / ((double)no_of_discrete_barriers)*(current_stock_price3 - initial_stock_price));
				mean_at_sampling_instant4[i] = initial_stock_price +
					(((double)i) / ((double)no_of_discrete_barriers)*(current_stock_price4 - initial_stock_price));
				variance_at_sampling_instant[i] = (((double)i) / ((double)no_of_discrete_barriers))
					*expiration_time*(1.0 - ((double)i) / ((double)no_of_discrete_barriers));
				probability_stock_path_has_not_hit_barrier_at_some_discrete_sampling_instant1 *=
					(1.0 - N((barrier_price - mean_at_sampling_instant1[i]) / sqrt(variance_at_sampling_instant[i])));
				probability_stock_path_has_not_hit_barrier_at_some_discrete_sampling_instant2 *=
					(1.0 - N((barrier_price - mean_at_sampling_instant2[i]) / sqrt(variance_at_sampling_instant[i])));
				probability_stock_path_has_not_hit_barrier_at_some_discrete_sampling_instant3 *=
					(1.0 - N((barrier_price - mean_at_sampling_instant3[i]) / sqrt(variance_at_sampling_instant[i])));
				probability_stock_path_has_not_hit_barrier_at_some_discrete_sampling_instant4 *=
					(1.0 - N((barrier_price - mean_at_sampling_instant4[i]) / sqrt(variance_at_sampling_instant[i])));
			}
			//calculate adjusted price by multipling probability_stock_path_has_not_hit_barrier_at_some_discrete_sampling_instant
			simulate_call_option_price_adjusted += (probability_stock_path_has_not_hit_barrier_at_some_discrete_sampling_instant1*breach_value_adjusted1*max(0.0, current_stock_price1 - strike_price) +
				probability_stock_path_has_not_hit_barrier_at_some_discrete_sampling_instant2*breach_value_adjusted2*max(0.0, current_stock_price2 - strike_price) +
				probability_stock_path_has_not_hit_barrier_at_some_discrete_sampling_instant3*breach_value_adjusted3*max(0.0, current_stock_price3 - strike_price) +
				probability_stock_path_has_not_hit_barrier_at_some_discrete_sampling_instant4*breach_value_adjusted4*max(0.0, current_stock_price4 - strike_price)) / 4.0;
			simulate_put_option_price_adjusted += (probability_stock_path_has_not_hit_barrier_at_some_discrete_sampling_instant1*breach_value_adjusted1*max(0.0, strike_price - current_stock_price1) +
				probability_stock_path_has_not_hit_barrier_at_some_discrete_sampling_instant2*breach_value_adjusted2*max(0.0, strike_price - current_stock_price2) +
				probability_stock_path_has_not_hit_barrier_at_some_discrete_sampling_instant3*breach_value_adjusted3*max(0.0, strike_price - current_stock_price3) +
				probability_stock_path_has_not_hit_barrier_at_some_discrete_sampling_instant4*breach_value_adjusted4*max(0.0, strike_price - current_stock_price4)) / 4.0;
	}
	//discount two prices to present
	simulate_call_option_price = exp(-risk_free_rate*expiration_time)*(simulate_call_option_price / ((double)no_of_trials));
	simulate_put_option_price = exp(-risk_free_rate*expiration_time)*(simulate_put_option_price / ((double)no_of_trials));
	simulate_call_option_price_adjusted = exp(-risk_free_rate*expiration_time)*(simulate_call_option_price_adjusted / ((double)no_of_trials));
	simulate_put_option_price_adjusted = exp(-risk_free_rate*expiration_time)*(simulate_put_option_price_adjusted / ((double)no_of_trials));
}


int main(int argc, char* argv[])
{
	std::chrono::time_point<std::chrono::system_clock> start, end;
	std::chrono::duration<double> elapsed_time;
	start = std::chrono::system_clock::now();//start calculating the time

	sscanf_s(argv[1], "%lf", &expiration_time);
	sscanf_s(argv[2], "%lf", &risk_free_rate);
	sscanf_s(argv[3], "%lf", &volatility);
	sscanf_s(argv[4], "%lf", &initial_stock_price);
	sscanf_s(argv[5], "%lf", &strike_price);
	sscanf_s(argv[6], "%d", &no_of_trials);
	sscanf_s(argv[7], "%d", &no_of_discrete_barriers);
	sscanf_s(argv[8], "%lf", &barrier_price);

	up_factor = exp(volatility*sqrt(expiration_time / ((double)no_of_discrete_barriers)));
	R = exp(risk_free_rate*expiration_time / ((double)no_of_discrete_barriers));
	uptick_prob = (R - (1 / up_factor)) / (up_factor - (1 / up_factor));

	cout << "European Down-and-Out Discrete Barrier Options Pricing via Monte Carlo Simulation" << endl;
	cout << "Expiration Time (Years) = " << expiration_time << endl;
	cout << "Risk Free Interest Rate = " << risk_free_rate << endl;
	cout << "Volatility (%age of stock value) = " << volatility * 100 << endl;
	cout << "Initial Stock Price = " << initial_stock_price << endl;
	cout << "Strike Price = " << strike_price << endl;
	cout << "Barrier Price = " << barrier_price << endl;
	cout << "Number of Trials = " << no_of_trials << endl;
	cout << "Number of Discrete Barriers = " << no_of_discrete_barriers << endl;
	cout << "--------------------------------------" << endl;
	cout << "Up Factor = " << up_factor << endl;
	cout << "Uptick Probability = " << uptick_prob << endl;
	cout << "--------------------------------------" << endl;
	adjusted_simulation();
	cout << "The average Call Price via explicit simulation of price paths             = " << simulate_call_option_price << endl;
	cout << "The average Call Price with Brownian-Bridge correction on the final price = " << simulate_call_option_price_adjusted << endl;
	cout << "--------------------------------------" << endl;
	cout << "The average Put Price via explicit simulation of price paths              = " << simulate_put_option_price << endl;
	cout << "The average Put Price with Brownian-Bridge correction on the final price  = " << simulate_put_option_price_adjusted << endl;
	cout << "--------------------------------------" << endl;

	end = std::chrono::system_clock::now();//stop calculating the time
	elapsed_time = end - start;//calculate the time used
	cout << "Elapsed time to complete: " << elapsed_time.count() << "s\n" << endl;//output the time
}