#include <iostream>
#include <fstream>
#include <vector>

using namespace std::vector;
using namespace std::cout;
using namespace std::cin;

//Calculate the residuals -- returns the form as a vector
vector<long double> CalcRes(vector<long double>& y_val, vector<long double>& model, int& size)
{
	long double res;
	vector<long double> res_vec;
	for(int i=0; i < size; ++i)
	{
		res = y_val(i)-model(i);
		res_vec.push_back(res);
	}
	return res_vec;
}

int main()
{
	// Choose the dataset

	// Choose the model

	// Ouput the stuff (where?)

	return 0;
}
