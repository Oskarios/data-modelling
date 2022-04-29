#include <iostream>
#include <math.h>
#include <fstream>
#include <vector>
#include <functional>

using std::vector;
using std::cout;
using std::cin;


long double MODEL_CexpMx(long double x, vector<long double> params, int num)
{
	long double result;
	cout << x << std::endl;
	result = params[0] * expl(params[1]/x);
	return result;
}



long double MODEL_POLY(long double& x, vector<long double>& params, int num)
{
	long double result = 0.0l;
	for(int i = 0; i < num; ++i)
	{
		result += params[i] * powl(x,i);
	}
	return result;
}

long double MODEL_LOG(long double x, vector<long double> params, int num)
{
	long double result;
	result = params[0] * log10(x) + params[1];
	return result;
}

long double LOG_MODEL_GRAD_dS_dm ()
{

}

long double LOG_MODEL_GRAD_dS_db ()
{

}

// Calculate d/dX (f) numerically -- return result in vector
// param: input function i.e. the residual, or the MSR
// param: differential variable
vector<long double> CalcDiff()
{

}


//Calculate the residuals -- returns the form as a vector
vector<long double> CalcRes(vector<long double>& y_val, vector<long double>& model, int& size)
{
	long double res;
	vector<long double> res_vec;
	for(int i=0; i < size; ++i)
	{
		res = y_val[i]-model[i];
		res_vec.push_back(res);
	}
	return res_vec;
}



long double CalcMSR(vector<long double>& residuals, int& size)
{
	long double MSR = 0.0l; // mean-squared-residual
	for(int i = 0; i < size; ++i)
	{
		// Technically the abs call is unnecessary bc r^2 > 0 no matter what
		// But I want this code to read/be better for sure
		MSR += abs(residuals[i]*residuals[i]);
	}
	MSR = MSR / size;
	return MSR;
}

// returns vector of y values predicted by model
vector<long double> MAP_MODEL(vector<long double>& x_vals, vector<long double> params, int set_size, int p_size, std::function<long double (long double, vector<long double>, int)> model)
{
	vector<long double> f_val;
	cout << "MAPPING \n";

	for (int i=0; i < set_size; ++i)
	{

		f_val.push_back(model(x_vals[i],params,p_size));
	}
	return f_val;
}


// Calculates gradient 2-vector for each parameter for exponential model
// Size refers to number of data, p_size refers to number of model parameters
// REQUIRES TESTING
vector<long double> EXP_MODEL_GRAD_dS_da(vector<long double>& residuals, vector<long double>& x_vals, vector<long double>& params, int size, int p_size)
{
	long double grad_j;
	vector<long double> grad;
	for(int j=0; j < p_size; ++j)
	{
		grad_j = 0.0l;
		for(int i=0; i < size; ++i)
		{
			grad_j += (-1.0l * MODEL_CexpMx(x_vals[i],params,2) * residuals[i])/params[j];
		}
		grad_j = 2.0l * grad_j / size;
		grad.push_back(grad_j);
	}
	return grad;
}

// Returns vector containing optimised model parameters (and final MSR) as calculated by steepest descent
// Start by implementing algorithm for just the small exponential model --> then expand for model choice
// init_param --> guesses for parameters, but should we have the lambda here? change later perhaps
// vector<long double> minimise_msr(vector<long double>& x_vals, vector<long double>& y_data, vector<long double>& init_param, int set_size)
// {
// 	// Let's assume that we're getting the number of parameters correct -> for each model it'll be 'hardcoded'
// 	int num_parameters = 2; // Recall that we're just starting with the exp model
// 	int max_iteration = 10000; // maximum number of iterations before giving up
//
// 	long double lambda = 0.001; // greed parameter
// 	long double delta_s; // change in S from step to step;
// 	// Calculate the vector?
//
// 	// Calculate the initial residuals, S, vector
// 	vector<long double> residuals = CalcRes()
//
// 	int i = 0;
// 	while(i < max_iteration) // we'll bring in the other criteria
// 	{
//
// 		++i;
// 	}
// }


//TEST
namespace DBG
{

	// Should make this more good later but good for now
	void DEBUG_PrintVec(vector<long double>& input_vec, int size)
	{
		for(int i = 0; i < size; ++i)
		{
			cout << "INPUT VEC; i=" << i << " VAL: " << input_vec[i] << std::endl;
		}
	}

	void DEBUG_PrintVec(vector<long double>& input_vec, vector<long double> expected, int size)
	{
		for(int i = 0; i < size; ++i)
		{
			cout << "INPUT VEC; i=" << i << " VAL: " << input_vec[i] << " EXP: "<< expected[i] << std::endl;
		}
	}

	void TEST_CalcRes(std::string test_name)
	{
		vector<long double> y_func = {1.1l, 1.9l, 3.1l};
		vector<long double> f_func = {1.0l, 2.0l, 3.0l};
		vector<long double> EXPECT = {0.1l, -0.1l, 0.1};
		int size = 3;

		vector<long double> TEST_RESULT = CalcRes(y_func, f_func, size);
		cout << std::endl << test_name << "\n\n";
		DEBUG_PrintVec(TEST_RESULT, EXPECT, size);
	}

	void TEST_MODEL_EXP(std::string test_name)
	{
		cout << "\n" << test_name << std::endl;

		long double expected, x;
		vector<long double> test_params;
		//TEST 1 C=2.5, m=1.5 (p0, p1)
		test_params={2.5l,1.5l};
		x = 1.1l;
		expected = 2.5l * expl(1.5/x);

		if(isnan(MODEL_CexpMx(x,test_params,2)))
		{
			cout << test_name << " ISNAN -- FAIL" << std::endl;
		}

		// cout << abs(MODEL_CexpMx(x,test_params,2) - expected) << std::endl;
		if(MODEL_CexpMx(x,test_params,2) == expected){
			cout << test_name << " 1 PASS" << std::endl;
		} else
		{
			cout << test_name << " 1 FAIL" << std::endl;
		}
		//TEST 2 -- we could require catching at x=0 and other similar cases

	}

	void TEST_MODEL_POLY(std::string test_name)
	{
		cout << "\n" << test_name << std::endl;
		vector<long double> t_pm;
		long double x;
		// TEST 1 : 4.2 + 2.1x + 3.4x^2 + 5x^3; x = 0.2
		t_pm = {4.2l,2.1l,3.4l,5.0l};
		x = 0.2l;
		long double expected = t_pm[0] + t_pm[1]*x + t_pm[2]*x*x + t_pm[3]*x*x*x;

		if(MODEL_POLY(x,t_pm,4) == expected)
		{
			cout << test_name << " 1 PASS" << std::endl;
		} else
		{
			cout << test_name << " 1 FAIL" << std::endl;
		}

		cout << MODEL_POLY(x,t_pm,4) << std::endl;
		cout << expected << std::endl;
	}

	void TEST_MODEL_LOG(std::string test_name)
	{
		cout << "\n" << test_name << std::endl;
		vector<long double> t_pm;
		long double x;
		//Test 1:
		t_pm = {3.13l,4.54l};
		x = 0.34l;

		long double expected = t_pm[0]*log10(x) + t_pm[1];
		long double MODEL = MODEL_LOG(x,t_pm,2);
		cout << MODEL << std::endl;
		cout << expected << std::endl;

		if(MODEL == expected)
		{
			cout << "PASS" << std::endl;
		}
	}

	void TEST_MAP_MODEL(std::string test_name)
	{
		cout << "\n" << test_name << std::endl;
		vector<long double> t_param = {1.2l,3.1l};
		vector<long double> x_vals = {1.0l,2.0l,3.0l,4.0l};
		// but what to do now?
		vector<long double> TEST_EXP = MAP_MODEL(x_vals, t_param, 4, 2, MODEL_CexpMx);
		vector<long double> TEST_LOG = MAP_MODEL(x_vals, t_param, 4, 2, MODEL_LOG);
		vector<long double> expon_expected;
		vector<long double> log_expected;

		for(int i = 0; i < 4; ++i)
		{
			expon_expected.push_back(MODEL_CexpMx(x_vals[i],t_param,2));
			log_expected.push_back(MODEL_LOG(x_vals[i],t_param,2));
		}

		cout << "TEST 1: Exponential Model" << std::endl;
		DEBUG_PrintVec(TEST_EXP,expon_expected, 4);
		cout << std::endl << "TEST 2: Log10 Model" << std::endl;
		DEBUG_PrintVec(TEST_LOG,log_expected, 4);
	}


	void RUN_TEST_HARNESS()
	{
		TEST_CalcRes("TEST: Calculating Residuals");

		TEST_MODEL_EXP("TEST: Calculating model");

		TEST_MODEL_POLY("TEST: Calculating polynomial model");

		TEST_MODEL_LOG("TEST: Calculating log10 model");

		TEST_MAP_MODEL("TEST: Map model function (f(x))");
	}
}



int main()
{
	DBG::RUN_TEST_HARNESS();


	// Choose the dataset

	// Choose the model

	// Ouput the stuff (where?)

	return 0;
}
