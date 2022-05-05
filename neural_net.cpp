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
	int index = 0;
	for(int i = 0; i < l_node_count[layer_index]; ++i)
	{
		for(int j = 0; j < l_node_count[layer_index-1]; ++j)
		{
			output_vec[i][j] = raw_weights[index];
			++index;
		}
	}
	return output_vec;
}

int main()
{

	vector<vector<double>> hidden_weights = {
		vector<double>(l_node_count[1]),
		vector<double>(l_node_count[0])
	};


	vector<vector<double>> output_weights = {
		vector<double>(l_node_count[2]),
		vector<double>(l_node_count[1])
	};

	// Set bias of all nodes to 0
	vector<double> hidden_bias(l_node_count[1],k_glob_bias);
	vector<double> output_bias(l_node_count[2],k_glob_bias); // -->

	vector<double> raw_hidden_weights = {0.4,0.2,0.1,-0.2,0.5,0.8};

	hidden_weights = set_weights(raw_hidden_weights,1,hidden_weights);

	// index = 0;
	vector<double> raw_output_weights = {0.2,-0.9,0.1};
	output_weights = set_weights(raw_output_weights,2,output_weights);

	vector<double> input = {1.0,2.0};

	vector<double> hidden_output(l_node_count[1]); // Calculate these by the weighting of the stuff

	double weight_i;
	double z;
	//Calculat the output of the sigmoid functions for the hidden layer
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

	std::cout << l_node_count.size() << std::endl;

	return 0;
}
