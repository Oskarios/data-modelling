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

double loss_function (vector<double>& output, vector<double>& actual)
{
	double loss = 0.0;
	int set_size = output.size();
	for(int i = 0; i < set_size; ++i)
	{
			loss += (output[i] - actual[i])*(output[i]-actual[i]);
	}
	loss = loss / set_size;
	return loss;
}

vector<vector<double>> prime_output_vec (vector<double>& input_set, vector<vector<double>>& output_vec)
{
	output_vec[0][0] = input_set[0];
	output_vec[0][1] = input_set[1];
	return output_vec;
}

// double network_result(vector<vector<double>>& neural_outputs)
// {
// 	double result = neural_outputs[2][0];
// 	return result;
// }

vector<double> grad_weights (vector<vector<double>>& weights, vector<double>& output_set, vector<double>& expected, vector<double>& input, vector<double>& grad_vec)
{
	double grad_j;
	// std::cout << "GRAD WEIGHTS "<< weights.size() << std::endl;
	for(int i = 0; i < weights.size(); ++i)
	{
		grad_j = 0.0;
		for(int j = 0; j < output_set.size(); ++j)
		{
			grad_j += (output_set[j] - expected[j]) * output_set[j]*(1-output_set[j]) * input[i];
		}
		grad_j = 2*grad_j / output_set.size();
		grad_vec[i] = grad_j;
	}
	return grad_vec;
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

	//Training set for neural network for OR gate
	vector<vector<double>> neural_inputs = {
		{0,0},
		{0,1},
		{1,0},
		{1,1}
	};

	vector<double> expected_or = {0,1,1,1};

	vector<double> neural_training_outputs(neural_inputs.size(),0.0);

	// vector<vector<vector<double>>> delta_weight = network;

	for(int i = 0; i < neural_inputs.size(); ++i)
	{
		layer_output = prime_output_vec(neural_inputs[i],layer_output);
		layer_output = feed_forward(network,net_bias,layer_output);
		neural_training_outputs[i] = layer_output[2][0];
		// std::cout << neural_inputs[i][0] << " "<< neural_inputs[i][1] << " " << neural_training_outputs[i] << "\n";
	}

	//Caculate the loss function
	double loss = loss_function(neural_training_outputs, expected_or);

	std::cout << "\nInitial Loss: " << loss << std::endl;

	// hoW the poop do I do this?
	double loss_new;
	double delta_loss = 0.01;

	vector<double> d_weight_hidden(l_node_count[1]*l_node_count[2],0);

	// std:: cout << "NETWORK " << network[1].size();
	std::cout << std::endl;

	// std::cout << "D_WEIGHT_HIDDEN  " << d_weight_hidden[0] << " "<< d_weight_hidden[1] <<  " " <<d_weight_hidden[2] << std::endl;

	while(true)
	{
		std::cout << "LOSS: " << loss << std::endl;
		// Set d_weight
		d_weight_hidden = grad_weights(network[0],neural_training_outputs,expected_or,layer_output[1],d_weight_hidden);
		// std::cout << "D_WEIGHT_HIDDEN  " << d_weight_hidden[0] << " "<< d_weight_hidden[1] <<  " " <<d_weight_hidden[2] << std::endl;

		//update the WEIGHTS
		for(int j = 0; j < d_weight_hidden.size(); ++j)
		{
			network[0][j][0] += 0.01 * d_weight_hidden[j];
		}

		// Calculate new output
		for(int i = 0; i < neural_inputs.size(); ++i)
		{
			layer_output = prime_output_vec(neural_inputs[i],layer_output);
			layer_output = feed_forward(network,net_bias,layer_output);
			neural_training_outputs[i] = layer_output[2][0];
			// std::cout << neural_inputs[i][0] << " "<< neural_inputs[i][1] << " " << neural_training_outputs[i] << "\n";
		}
		// Calculate new loss
		loss_new = loss_function(neural_training_outputs, expected_or);
		//repeat?
		delta_loss = loss_new - loss;
		loss = loss_new;
	}


	// std::cout << layer_output[2][0] << std::endl;


	return 0;
}
