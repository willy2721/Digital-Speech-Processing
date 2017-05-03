#include "hmm.h"
#include <math.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <iterator>
#include <map>
#include <algorithm>
using namespace std;

namespace patch
{
    template < typename T > std::string to_string( const T& n )
    {
        std::ostringstream stm ;
        stm << n ;
        return stm.str() ;
    }
}


int main(int argc, char *argv[])
{
	if(argc != 5){
		printf("%s unable to execute due to wrong number of arguments", argv[0]);
	}
	
	int rounds = atoi(argv[1]);
	char* init_model = argv[2];
	char* sequence_model = argv[3];
	char* out_model = argv[4];
	

	// Set up global variables//
	int time_period = 0;
	int state_num = 0;
	int observ_num = 0;
	int sample_num = 0;

	// Define hmm_initial.transition map to map chars to int //
	map<char,int> seq_map;
	seq_map['A'] = 0;
	seq_map['B'] = 1;
	seq_map['C'] = 2;
	seq_map['D'] = 3;
	seq_map['E'] = 4;
	seq_map['F'] = 5;

	// Load seq_model //
	vector<string> seq_model;
	string line;
  	ifstream myfile (sequence_model);
  	if (myfile.is_open()){
    	while ( getline (myfile,line) ){
    		if(time_period == 0) time_period = line.length();
    		seq_model.push_back(line);
    		sample_num++;
    	}
    	myfile.close();
  	}

  	// Create integer representation of seq_model (A = 0 ... F = 5) // 
  	vector< vector<int> > seq_int;
  	for(int i = 0; i < seq_model.size(); i++){
  		vector<int> tmp;
  		for(int j = 0; j < seq_model[i].size(); j++){
  			tmp.push_back(seq_map.find(seq_model[i][j])->second);
  		}
  		seq_int.push_back(tmp);
  	}

  	// Load initial HMM //
	HMM hmm_initial;
	loadHMM( &hmm_initial, init_model);
	
	// Save initial HMM as pi, hmm_initial.transition and hmm_initial.observation;
	state_num = hmm_initial.state_num;
	observ_num = hmm_initial.observ_num;


	/****************************************************************************/

	for(int iter = 0; iter < rounds; iter++){
		
		// Create variables to store accumulation
	  	// gamma_one to store the sum of gamma ones for the hmm_initial.initial update
		// gamma_t_min to store the sums of denominators for the a_ij update
		// gamma_t to store the sums of denominators for the b_jk update
		vector<double> gamma_one(state_num,0);
		vector<double> gamma_t_min(state_num,0);
		vector<double> gamma_t(state_num,0);

		// epsilon_t_min to store the sum of nominators for the a_ij update
		// gamma_t_each to store the sum of nominators for the b_jk update
		vector<vector<double> > epsilon_t_min(state_num,vector<double>(state_num,0));
		vector<vector<double> > gamma_t_each(observ_num,vector<double>(state_num,0));


		// Repeat for number of samples
		for(int n = 0; n < sample_num; n++){
			// Create 2D vector alpha, beta and gamma //
		  	vector< vector<double> > alpha, beta, gamma;
		  	vector<double> beta_value, alpha_value, gamma_value;

		  	// Calculate the initial values
		  	for(int i = 0; i < state_num; i++){
		  		alpha_value.push_back(hmm_initial.initial[i] * hmm_initial.observation[i][seq_int[n][0]]); // FIX!!  LEFT 0 -> ONLY FIRST LINE
		  		beta_value.push_back(1);
		  	}
		  	alpha.push_back(alpha_value);
		  	alpha_value.clear();
		  	beta.push_back(beta_value);
		  	beta_value.clear();

		  	// Calculate values for each time period after initial alpha (forward) and initial beta (backward)
		  	for(int t = 1; t < time_period; t++){
		  		// For each alpha[i][j]
		  		for(int i = 0; i < state_num; i++){	
		  			double tmp_sum_alpha = 0;
		  			double tmp_b = 0;
		  			double tmp_sum_beta = 0;
		  			for(int j = 0; j < state_num; j++){
		  				tmp_sum_alpha += alpha[t - 1][j] * hmm_initial.transition[j][i];
		  				tmp_sum_beta += hmm_initial.transition[i][j] * hmm_initial.observation[seq_int[n][time_period - t]][j] * beta[t-1][j]; // FIX!!
		  			}
		  			beta_value.push_back(tmp_sum_beta);
		  			alpha_value.push_back(tmp_sum_alpha * hmm_initial.observation[seq_int[n][t]][i]); // FIX!!	
		  		}
		  		// Store the rows for alpha and beta
		  		alpha.push_back(alpha_value);
		  		alpha_value.clear();
		  		beta.push_back(beta_value);
		  		beta_value.clear();
		  	}
		  	reverse(beta.begin(),beta.end());


			// Calculate values for each time period for gamma 
		  	for(int t = 0; t < time_period; t++){
		  		double denom = 0;
		  		for(int i = 0; i < state_num; i++){
		  			denom += alpha[t][i] * beta[t][i];
		  		}
		  		for(int i = 0; i < state_num; i++){
		  			gamma_value.push_back(alpha[t][i] * beta[t][i] / denom);
		  		}
		  		gamma.push_back(gamma_value);
		  		gamma_value.clear();
		  	}

			// Create 3D vector epsilon
		  	vector<vector<vector<double> > > epsilon;
		  	vector<vector<double> > epsilon_ij;
		  	vector<double> epsilon_value;

		  	// Calculate values for each time period for epsilon
		  	for(int t = 0; t < time_period - 1; t++){
		  		double denom = 0;
		  		for(int i = 0; i < state_num; i++){
		  			for(int j = 0; j < state_num; j++){
		  				denom += alpha[t][i] * hmm_initial.transition[i][j] * hmm_initial.observation[seq_int[n][t + 1]][j] * beta[t + 1][j]; // FIX!!
		  			}	
		  		}
		  		for(int i = 0; i < state_num; i++){
		  			for(int j = 0; j < state_num; j++){
		  				epsilon_value.push_back(alpha[t][i] * hmm_initial.transition[i][j] * hmm_initial.observation[seq_int[n][t + 1]][j] * beta[t + 1][j] / denom); // FIX!!
		  			}
		  			epsilon_ij.push_back(epsilon_value);
		  			epsilon_value.clear();
		  		}
		  		epsilon.push_back(epsilon_ij);
		  		epsilon_ij.clear();
		  	}

		  	// Accumulate gamma_one, epsilon_t_min, gamma_t_min, gamma_tk, gamma_t ...
		  	for(int i = 0; i < state_num; i++){
		  		gamma_one[i] += gamma[0][i];
		  	}

		  	for(int t = 0; t < time_period; t++){
		  		// for epsilon_t_min and gamma_t_min
		  		if(t != time_period - 1){
		  			for(int i = 0; i < state_num; i++){
		  				gamma_t_min[i] += gamma[t][i];
		  				gamma_t[i] += gamma[t][i];
		  				gamma_t_each[seq_int[n][t]][i] += gamma[t][i]; // FIX!!
		  				for(int j = 0; j < state_num; j++){
		  					epsilon_t_min[i][j] += epsilon[t][i][j];
		  				}
		  			}
		  		}
		  		else{
				// deal with last time period (gamma_t_each and gamma_t)
		  			for(int i = 0; i < state_num; i++){
		  				gamma_t[i] += gamma[t][i];
		  				gamma_t_each[seq_int[n][t]][i] += gamma[t][i]; // FIX!!
		  			}  			
		  		}
		  	}
		}


	  	// Update hmm_initial.initial, hmm_initial.transition and hmm_initial.observation
	  	for(int i = 0; i < state_num; i++){
	  		hmm_initial.initial[i] = gamma_one[i] / sample_num; // FIX!!
	  		for(int j = 0; j < state_num; j++){
	  			hmm_initial.transition[i][j] = epsilon_t_min[i][j] / gamma_t_min[i];
	  		}
	  		for(int k = 0; k < observ_num; k++){
	  			hmm_initial.observation[k][i] = gamma_t_each[k][i] / gamma_t[i];
	  		}
	  	}
	  	
	}

	// Write hmm_initial.initial, hmm_initial.transition, hmm_initial.observation to file
	ofstream outfile;
	outfile.open(argv[4]);
	string output = "";
	output += ("initial: " + string(patch::to_string(state_num)) + "\n");
	for(int i = 0; i < state_num; i++){
		output += patch::to_string(hmm_initial.initial[i]);
		if(i < state_num - 1)
			output += " ";
	}
	output += "\n\n";
	output += ("transition: " + string(patch::to_string(state_num)) + "\n");
	for(int i = 0; i < state_num; i++){
		for(int j = 0; j < state_num; j++){
			output += patch::to_string(hmm_initial.transition[i][j]);
			if(j < state_num - 1)
				output += " ";
		}
		output += "\n";
	}
	output += "\n";
	output += ("observation: " + string(patch::to_string(observ_num)) + "\n");
	for(int i = 0; i < observ_num; i++){
		for(int j = 0; j < state_num; j++){
			output += patch::to_string(hmm_initial.observation[i][j]);
			if(j < state_num - 1)
				output += " ";
		}
		if(i < observ_num - 1)
			output += "\n";
	}

	outfile << output << endl;
	outfile.close();
		

  	/* Test print 
  	// Printing hmm_initial.initial, hmm_initial.transition and hmm_initial.observation
	printf("hmm_initial.initial:\n");
  	for(int i = 0; i < state_num; i++){
  		printf("%f ", hmm_initial.initial[i]);
  	}
  	printf("\n");

  	printf("hmm_initial.transition:\n");
  	for(int i = 0; i < state_num; i++){
  		for(int j = 0; j < state_num; j++){
  			printf("%f ", hmm_initial.transition[i][j]);	
  		}
  		printf("\n");
  	}
  	printf("\n");

  	printf("hmm_initial.observation:\n");
  	for(int i = 0; i < observ_num; i++){
  		for(int j = 0; j < state_num; j++){
  			printf("%f ", hmm_initial.observation[i][j]);	
  		}
  		printf("\n");
  	}
  	printf("\n");

  	printf("\n");
  	*/

  	/* Test print
  	// Printing accumulations
	printf("gamma_one:\n");
  	for(int i = 0; i < state_num; i++){
	  	printf("%f ",gamma_one[i]);
	}
	printf("\n\n");

  	printf("epsilon_t_min:\n");
  	for(int i = 0; i < state_num; i++){
	  	for(int j = 0; j < state_num; j++){
	  		printf("%f ",epsilon_t_min[i][j]);
	  	}
	  	printf("\n");
	}
	printf("\n");

	printf("gamma_t_min:\n");
  	for(int i = 0; i < state_num; i++){
	  	printf("%f ",gamma_t_min[i]);
	}
	printf("\n\n");

	printf("gamma_t:\n");
  	for(int i = 0; i < state_num; i++){
	  	printf("%f ",gamma_t[i]);
	}
	printf("\n\n");

	printf("gamma_t_each:\n")
	for(int i = 0; i < observ_num; i++){
		printf("gamma %i:\n", i);
		for(int j = 0; j < state_num; j++){
			printf("%f ",gamma_t_each[i][j]);	
		}	
		printf("\n");		
	}
  	
  	/* Test print
	// Number of time periods
  	printf("%d\n",time_period);
	// Initial alpha
	for(int i = 0; i < hmm_initial.state_num; i++){
  		printf("%f ", alpha[0][i]);
  	}
  	// Full alpha (or beta)
  	for(int i = 0; i < state_num; i++){
  		for(int j = 0; j < state_num; j++){
  			printf("%f ",alpha[i][j]);
  		}
  		printf("\n");
  	}
	// Sequence in integers
  	for(int i = 0; i < 10; i++){
  		for(int j = 0; j < 10; j++){
  			printf("%d ", seq_int[i][j]);
  		}
  		printf("\n");
  	}
  	// Epsilon
  	for(int t = 0; t < 1; t++){
  		double sum = 0;
	  	for(int i = 0; i < state_num; i++){
	  		for(int j = 0; j < state_num; j++){
	  			sum += epsilon[t][i][j];
	  			printf("%f ",epsilon[t][i][j]);
	  		}
	  		printf("\n");
	  	}
	  	printf("%f\n",sum);
	}
	*/
		

	return 0;
}
