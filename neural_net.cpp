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

		// std::cout << "VTEMP: ";
		// for(int k = 0; k < vtemp.size(); ++k)
		// {
		// 	std::cout << vtemp[k] << "\t";
		// }
		// std::cout << std::endl;

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
	// std::cout << "Init test: " << hidden_weights[0][0] << std::endl;

	// Set bias of all nodes to 0
	vector<double> hidden_bias(l_node_count[1],0);
	vector<double> output_bias(l_node_count[2],0); // -->

	vector<double> raw_hidden_weights = {0.4,0.2,0.1,-0.2,0.5,0.8};

	hidden_weights = set_weights(raw_hidden_weights,1,hidden_weights);

	// std::cout << "We reached here\n";

	// for(int i = 0; i < 3; ++i)
	// {
	// 	for(int j = 0; j < 2; ++j)
	// 	{
	// 		std::cout << i << "\t" << j << "\t" << hidden_weights[i][j] << "\n";
	// 	}
	// }
	//

	// index = 0;
	vector<double> raw_output_weights = {0.2,-0.9,0.1};
	output_weights = set_weights(raw_output_weights,2,output_weights);

	vector<vector<double>> layer_output;

	vector<double> n_input = {1.0,2.0};
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

	for(int i = 1; i < layers; ++i)
	{
		for(int j = 0; j < l_node_count[i]; ++j)
		{
			weight_i = 0.0;
			for(int k = 0; k < l_node_count[i-1]; ++k)
			{
				// std::cout << i << j << k << std::endl;
				weight_i += network[i-1][j][k] * layer_output[i - 1][k];
			}

			// std::cout << "BIAS: " << net_bias[i-1][j] << "\n";
			z = net_bias[i-1][j] + weight_i;
			layer_output[i][j] = activation(z);
			std::cout << "OUTPUT: " << layer_output[i][j] << "\n";

			// layer_output[i+1].push_back(activation(z));

			// std::cout << i << "\t" << layer_output[i][j] << "\n";
		}
	}

	std::cout << layer_output[2][0] << std::endl;


	return 0;
}
