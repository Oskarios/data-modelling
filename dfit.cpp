#include <iostream>
#include <math.h>
#include <fstream>
#include <vector>
#include <functional>
#include <iterator>
#include <sstream>
#include <string>

using std::vector;
using std::cout;
using std::cin;

// enum Model {};
const double k_GRAD_VEC_TOL = 0.000001;
const double k_DELTA_MSR_TOL = 0.000001;

// Returns the line from the geometry file as a vector
vector<double> SplitLine(std::string line){
    std::istringstream iss(line);

    return vector< double>{
        std::istream_iterator<double>(iss),
        std::istream_iterator<double>()
    };
}


double MODEL_CexpMx(double x, vector<double> params, int num)
{
	double result;
	// cout << x << std::endl;
	result = params[0] * exp(params[1]/x);
	return result;
}



double MODEL_POLY(double& x, vector<double>& params, int num)
{
	double result = 0.0l;
	for(int i = 0; i < num; ++i)
	{
		result += params[i] * pow(x,i);
	}
	return result;
}

double MODEL_LOG(double x, vector<double> params, int num)
{
	double result;
	result = params[0] * log10(x) + params[1];
	return result;
}

double LOG_MODEL_GRAD_dS_dm ()
{

}

double LOG_MODEL_GRAD_dS_db ()
{

}

// Calculate d/dX (f) numerically -- return result in vector
// param: input function i.e. the residual, or the MSR
// param: differential variable
vector<double> CalcDiff()
{

}


//Calculate the residuals -- returns the form as a vector
vector<double> CalcRes(vector<double>& y_val, vector<double>& model, int& size)
{
	double res;
	vector<double> res_vec;
	for(int i=0; i < size; ++i)
	{
		res = y_val[i]-model[i];
		res_vec.push_back(res);
	}
	return res_vec;
}



double CalcMSR(vector<double>& residuals, int& size)
{
	double MSR = 0.0; // mean-squared-residual
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
vector<double> MAP_MODEL(

		vector<double>& x_vals,
		vector<double>& params, 
		int set_size, 
		int p_size, 
		std::function<double (double, vector<double>, int)> model
	
	)
{
	vector<double> f_val;
	// cout << "MAPPING \n";
	for (int i=0; i < set_size; ++i)
	{
		f_val.push_back(model(x_vals[i],params,p_size));
	}
	return f_val;
}


// Calculates gradient 2-vector for each parameter for exponential model
// Size refers to number of data, p_size refers to number of model parameters
// REQUIRES TESTING
vector<double> EXP_MODEL_GRAD_dS_da(vector<double>& residuals, vector<double>& x_vals, vector<double>& params, int size, int p_size)
{
	double grad_j;
	vector<double> grad;
	vector<double> sum_term_param; // To contain {C,x_i}
	sum_term_param.push_back(params[0]); // contains {C} now
	for(int j=0; j < p_size; ++j)
	{
		grad_j = 0.0;
		for(int i=0; i < size; ++i)
		{	
			sum_term_param.push_back(x_vals[i]);
			grad_j += (-1.0 * MODEL_CexpMx(x_vals[i],params,2) * residuals[i])/sum_term_param[j]; //don't use this bit
			sum_term_param.pop_back();
		}
		grad_j = 2.0 * grad_j / size;
		grad.push_back(grad_j);
	}
	return grad;
}

vector<double> update_parameters (vector<double>& parameters, vector<double>& grad, int& dimension, double& greed)
{
	vector<double> new_params;
	for(int i=0; i  < dimension; ++i)
	{
		new_params.push_back(parameters[i] - grad[i] * greed);
	}
	return new_params;
}

double grad_mod_squared(vector<double>& grad_vector, int dimension)
{
	double gms = 0.0;
	for(int i = 0; i < dimension; ++i)
	{
		gms += abs(grad_vector[i])*abs(grad_vector[i]);
	}
	return gms;
}


// Returns vector containing optimised model parameters (and final MSR) as calculated by steepest descent
// Start by implementing algorithm for just the small exponential model --> then expand for model choice
// init_param --> guesses for parameters, but should we have the lambda here? change later perhaps
vector<double> minimise_msr(
	vector<double>& x_vals, 
	vector<double>& y_data, 
	vector<double>& params, 
	int set_size, 
	int num_parameters,
	std::function<double (double, vector<double>, int)> model_fn,
	std::function<vector<double> (vector<double>&, vector<double>&, vector<double>&, int, int)> grad_vec_fn,
	std::ofstream& output_file
	)
{

	int max_iteration = 100000;// maximum number of iterations before giving up


	double lambda = 0.001; // greed parameter
	double delta_s = 1000; // change in S from step to step;
	double gms = 1000; 	  // initial loop values	
	


	// Calculate the initial residuals, S, vector
	vector<double> f_vals = MAP_MODEL(x_vals,params,set_size,num_parameters,model_fn);
	vector<double> residuals = CalcRes(y_data,f_vals,set_size);

	double MSR_old = CalcMSR(residuals, set_size);
	// Calculate the vector?

	vector<double> grad_vector = grad_vec_fn(residuals,x_vals,params,set_size,num_parameters);
	double MSR_new;

	int i = 0;
	while(!(delta_s < k_DELTA_MSR_TOL && gms < k_GRAD_VEC_TOL)) // we'll bring in the other criteria
	{
		// Change the parameters based on the gradient vector
		params = update_parameters(params, grad_vector, num_parameters, lambda);
		// Calculate the model with the updated model parameters

		f_vals = MAP_MODEL(x_vals,params,set_size,num_parameters,model_fn);
		residuals.clear();
		residuals = CalcRes(y_data,f_vals,set_size);
		MSR_new = CalcMSR(residuals,set_size);

		delta_s = MSR_new - MSR_old;
		MSR_old = MSR_new;
		grad_vector.clear();
		grad_vector = grad_vec_fn(residuals,x_vals,params,set_size,num_parameters);

		gms = grad_mod_squared(grad_vector,num_parameters);
		output_file << params[0] << "	" << params[1] << "	 " << delta_s << "	" << gms << "	" << MSR_new << std::endl;

		++i;
		if(i >= max_iteration){
			break;
		}
	}

	cout << "\nIterations counted: " << i << std::endl;

	vector<double> minimised;
	for(int i = 0; i < num_parameters; ++i)
	{
		minimised.push_back(params[i]);
	}
	minimised.push_back(MSR_new);
	return minimised;

}


//TEST
namespace DBG
{

	// Should make this more good later but good for now
	void DEBUG_PrintVec(vector<double>& input_vec, int size)
	{
		for(int i = 0; i < size; ++i)
		{
			cout << "INPUT VEC; i=" << i << " VAL: " << input_vec[i] << std::endl;
		}
	}

	void DEBUG_PrintVec(vector<double>& input_vec, vector<double> expected, int size)
	{
		for(int i = 0; i < size; ++i)
		{
			cout << "INPUT VEC; i=" << i << " VAL: " << input_vec[i] << " EXP: "<< expected[i] << std::endl;
		}
	}

	void TEST_CalcRes(std::string test_name)
	{
		vector<double> y_func = {1.1, 1.9, 3.1};
		vector<double> f_func = {1.0, 2.0, 3.0};
		vector<double> EXPECT = {0.1, -0.1, 0.1};
		int size = 3;

		vector<double> TEST_RESULT = CalcRes(y_func, f_func, size);
		cout << std::endl << test_name << "\n\n";
		DEBUG_PrintVec(TEST_RESULT, EXPECT, size);
	}

	void TEST_MODEL_EXP(std::string test_name)
	{
		cout << "\n" << test_name << std::endl;

		double expected, x;
		vector<double> test_params;
		//TEST 1 C=2.5, m=1.5 (p0, p1)
		test_params={2.5,1.5};
		x = 1.1;
		expected = 2.5 * expl(1.5/x);

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
		vector<double> t_pm;
		double x;
		// TEST 1 : 4.2 + 2.1x + 3.4x^2 + 5x^3; x = 0.2
		t_pm = {4.2,2.1,3.4,5.0};
		x = 0.2l;
		double expected = t_pm[0] + t_pm[1]*x + t_pm[2]*x*x + t_pm[3]*x*x*x;

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
		vector<double> t_pm;
		double x;
		//Test 1:
		t_pm = {3.13,4.54};
		x = 0.34;

		double expected = t_pm[0]*log10(x) + t_pm[1];
		double MODEL = MODEL_LOG(x,t_pm,2);
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
		vector<double> t_param = {1.0,3.1};
		vector<double> x_vals = {1.0,2.0,3.0,4.0};
		// but what to do now?
		vector<double> TEST_EXP = MAP_MODEL(x_vals, t_param, 4, 2, MODEL_CexpMx);
		vector<double> TEST_LOG = MAP_MODEL(x_vals, t_param, 4, 2, MODEL_LOG);
		vector<double> expon_expected;
		vector<double> log_expected;

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

	
	const int k_num_models = 1;

	const vector<std::function<double (double, vector<double>, int)>> k_fn_arr = {MODEL_CexpMx,MODEL_LOG};
	const vector<std::function<vector<double> (vector<double>&, vector<double>&, vector<double>&, int, int)>> k_grad_fn_arr = {EXP_MODEL_GRAD_dS_da};
	vector<vector<double>> k_param_init = {{10.0,10.0}}; // would like this to be a constant but it doesn't like it --> I am not smart enough for pointers etc.
	const vector<int> k_model_param_dim = {2,2};
	const vector<std::string> k_in_fnames = {"data_in/simple_set.dat"};
	const vector<std::string> k_out_fname = {"data_out/simple_set_out.dat"};


	//iterate over the datasets/models
	std::ifstream ifile;
	std::ofstream ofile;

	vector<double> x_vals;
	vector<double> y_vals;
	vector<double> minimised_params;
	int lines;

	for(int i = 0; i < k_num_models; ++i)
	{
		minimised_params.clear();
		//Open the relevant file
		ifile.open(k_in_fnames[i]);
		ofile.open(k_out_fname[i],std::ios::trunc);

		x_vals.clear();
		y_vals.clear();

		lines = 0;

		for(std::string line; getline(ifile,line);)
		{
			if(line.front() == '#')
			{
				cout << "Line skipped " << line << std::endl;
				continue;
			}

			x_vals.push_back(SplitLine(line)[0]);
			y_vals.push_back(SplitLine(line)[1]);
			++lines;
		}

		ifile.close();

		minimised_params = minimise_msr(x_vals,y_vals,k_param_init[i],lines,k_model_param_dim[i],k_fn_arr[i],k_grad_fn_arr[i],ofile);
		ofile.close();

		cout << std::endl;
		for(int l = 0; l < k_model_param_dim[i]; ++l)
		{
			cout << "P" << l << "	";
		}
		cout << "MSR" << std::endl;
		for(int k = 0; k <= k_model_param_dim[i]; ++k)
		{
			cout << minimised_params[k] << "	";
		}

		cout << std::endl;

	}

	return 0;
}
