#include <iostream>
#include <math.h>
#include <fstream>
#include <vector>

using std::vector;
using std::cout;
using std::cin;


long double MODEL_CexpMx(long double& x, vector<long double>& params, int num)
{
	long double result;
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
		//TEST 2

	}

	void TEST_MODEL_POLY(std::string test_name)
	{

	}

}



int main()
{
	DBG::TEST_CalcRes("TEST: Calculating Residuals");

	DBG::TEST_MODEL_EXP("TEST: Calculating model");

	// Choose the dataset

	// Choose the model

	// Ouput the stuff (where?)

	return 0;
}
