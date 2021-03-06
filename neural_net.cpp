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
			++index;
		}
		output_vec[i].clear();
		output_vec[i] = vtemp;
	}
	return output_vec;
}

vector<vector<double>> feed_forward (
	vector<vector<vector<double>>>& weights,
	vector<vector<double>>& node_bias,
	vector<vector<double>>& node_outputs
  )
{
	double weight_i;
	double z;
	for(int i = 1; i < layers; ++i)
	{
		for(int j = 0; j < l_node_count[i]; ++j)
		{
			weight_i = 0.0;
			for(int k = 0; k < l_node_count[i-1]; ++k)
			{
				weight_i += weights[i-1][j][k] * node_outputs[i-1][k];
			}
			z = node_bias[i-1][j] + weight_i;
			node_outputs[i][j] = activation(z);
		}
	}

	return node_outputs;
}

int main()
{

	vector<vector<double>> hidden_weights(3,vector<double>(2,0)); // size 2
	std::cout << hidden_weights.size() << "\n";
	vector<vector<double>> output_weights(1,vector<double>(3,0));
	// std::cout << "Init test: " << hidden_weights[0][0] << std::endl;

	// Set bias of all nodes to 0
	vector<double> hidden_bias(l_node_count[1],0);
	vector<double> output_bias(l_node_count[2],0); // -->

	vector<double> raw_hidden_weights = {0.4,0.2,0.1,-0.2,0.5,0.8};

	hidden_weights = set_weights(raw_hidden_weights,1,hidden_weights);

	// index = 0;
	vector<double> raw_output_weights = {0.2,-0.9,0.1};
	output_weights = set_weights(raw_output_weights,2,output_weights);

	vector<vector<double>> layer_output;

	vector<double> n_input = {2.0,1.0};
	vector<double> hidden_output(l_node_count[1],0.0); // Calculate these by the weighting of the stuff
	// vector<double> n_output;
	vector<double> output_output(l_node_count[2],0.0);

	layer_output.push_back(n_input); 	   // i = 0
	layer_output.push_back(hidden_output); // i = 1
	layer_output.push_back(output_output); // i = 2

	vector<vector<vector<double>>> network;
	network.push_back(hidden_weights);
	network.push_back(output_weights);

	vector<vector<double>> net_bias;
	net_bias.push_back(hidden_bias);
	net_bias.push_back(output_bias);

	double weight_i;
	double z;
	//Calculat the output of the sigmoid functions for the hidden layer

	// std::cout << layer_output[0][0] << "\t" << layer_output[0][1] << "\n";
	// std::cout << layer_output[1][0] << "\t" << layer_output[2][3] << "\n";

	// Iterate over the layers of weights of the network etc.
	// std::cout << "RAH" << network[1][1][0];
	std::cout << "\n";

	layer_output = feed_forward(network,net_bias,layer_output);


	std::cout << "NN INPUT: " << layer_output[0][0] << "\t" << layer_output[0][1] << "\n";

	std::cout << "NN RESULT: ";
	std::cout << layer_output[2][0] << std::endl;


	return 0;
}
