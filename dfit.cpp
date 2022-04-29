#include <iostream>
#include <fstream>
#include <vector>

using std::vector;
using std::cout;
using std::cin;

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

int main()
{
	TEST_CalcRes("TEST: Calculating Residuals");

	// Choose the dataset

	// Choose the model

	// Ouput the stuff (where?)

	return 0;
}
