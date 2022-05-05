#include <iostream>
#include <vector>
#include <cmath>

using std::vector;

const int layers = 3;
const vector<int> l_node_count = {2,3,1}; // Number of nodes per layer
const double k_glob_bias = 0.0;


double activation (double z)
{
	return 1.0/(1.0+exp(-1.0*z));
}

vector<vector<double>> set_weights(vector<double>& raw_weights, int layer_index, vector<vector<double>>& output_vec)
{
	vector<double> vtemp;

	int index = 0;
	for(int i = 0; i < l_node_count[layer_index]; ++i)
	{
		vtemp.clear();
		for(int j = 0; j < output_vec[i].size(); ++j)
		{
			vtemp.push_back(raw_weights[index]);
			// std::cout << raw_weights[index] << "\t" << i << "\t" << j << "\t" << output_vec[j][i] << std::endl;
			++index;
		}

		std::cout << "VTEMP: ";
		for(int k = 0; k < vtemp.size(); ++k)
		{
			std::cout << vtemp[k] << "\t";
		}
		std::cout << std::endl;

		output_vec[i].clear();
		output_vec[i] = vtemp;
	}
	return output_vec;
}

int main()
{

	vector<vector<double>> hidden_weights(3,vector<double>(2,0)); // size 2
	vector<double> random;
	std::cout << hidden_weights.size() << "\n";


	vector<vector<double>> output_weights(1,vector<double>(3,0));


	std::cout << "Init test: " << hidden_weights[0][0] << std::endl;

	// Set bias of all nodes to 0
	vector<double> hidden_bias(l_node_count[1],0);
	vector<double> output_bias(l_node_count[2],0); // -->

	vector<double> raw_hidden_weights = {0.4,0.2,0.1,-0.2,0.5,0.8};

	hidden_weights = set_weights(raw_hidden_weights,1,hidden_weights);

	std::cout << "We reached here\n";

	for(int i = 0; i < 3; ++i)
	{
		for(int j = 0; j < 2; ++j)
		{
			std::cout << i << "\t" << j << "\t" << hidden_weights[i][j] << "\n";
		}
	}


	// index = 0;
	vector<double> raw_output_weights = {0.2,-0.9,0.1};
	output_weights = set_weights(raw_output_weights,2,output_weights);

	vector<double> input = {1.0,2.0};

	vector<double> hidden_output; // Calculate these by the weighting of the stuff


	double weight_i;
	double z;
	//Calculat the output of the sigmoid functions for the hidden layer

	std::cout << "Output from hidden layer: \n";
	for(int i = 0; i < l_node_count[1]; ++i)
	{
		weight_i = 0.0;
		for(int j = 0; j < l_node_count[1-1]; ++j)
		{
			// Sum up your weights
			weight_i += hidden_weights[i][j] * input[j];
		}
		// Add the bias for each node
		z = hidden_bias[i] + weight_i;
		// Calculate the sigmoid/switching function
		hidden_output.push_back(activation(z));
		// Assign the hidden node output for node i (in hidden layer)
		//Debug line
		std::cout << i << "\t" << hidden_output[i] << std::endl;
	}

	// std::cout << l_node_count.size() << std::endl;

	vector<double> output_output(l_node_count[2]);
	std::cout << "\nNetwork Output:\n";

	// PROCESS THE OUTPUT LAYER
	for(int i = 0; i < l_node_count[2]; ++i)
	{
		weight_i = 0.0;
		for(int j = 0; j < l_node_count[1]; ++j)
		{
			// std::cout << output_weights[i][j] << "\t" << hidden_output[j] << std::endl;

			weight_i += output_weights[i][j] * hidden_output[j];
		}

		z = output_bias[i] + weight_i;

		output_output.clear();
		output_output.push_back(activation(z));
		std::cout << i << "\t" << output_output[i] << std::endl;
	}



	return 0;
}
