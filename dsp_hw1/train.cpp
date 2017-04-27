#include "hmm.h"
#include <math.h>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <iterator>
#include <map>
#include <algorithm>
using namespace std;

int main()
{
/*
	HMM hmms[5];
	load_models( "modellist.txt", hmms, 5);
	dump_models( hmms, 5);
*/

	// Set up global variables//
	int time_period = 0;
	vector<double> pi;
	vector<vector<double> > a, b;
	int state_num = 0;
	int observ_num = 0;

	// Define a map to map chars to int //
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
  	ifstream myfile ("seq_model_01.txt");
  	if (myfile.is_open()){
    	while ( getline (myfile,line) ){
    		if(time_period == 0) time_period = line.length();
    		seq_model.push_back(line);
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
	loadHMM( &hmm_initial, "model_init.txt" );
	
	// Save initial HMM as pi, a and b;
	state_num = hmm_initial.state_num;
	observ_num = hmm_initial.observ_num;

	// Save hmm_initial.initial as pi
	for(int i = 0 ; i < state_num; i++){
		pi.push_back(hmm_initial.initial[i]);
	}

	// Save hmm_initial.transition as a
	for(int i = 0; i < state_num; i++){
		vector<double> tmp;
		for(int j = 0; j < state_num; j++){
			tmp.push_back(hmm_initial.transition[i][j]);
		}
		a.push_back(tmp);
		tmp.clear();	
	}

	// Save hmm_initial.observation as b
	for(int i = 0; i < observ_num; i++){
		vector<double> tmp;
		for(int j = 0; j < state_num; j++){
			tmp.push_back(hmm_initial.observation[i][j]);
		}
		b.push_back(tmp);
		tmp.clear();
	}

  	// Create 2D vector alpha //
  	vector< vector<double> > alpha;
  	vector<double> alpha_init, alpha_value;

  	// Calculate the initial values
  	for(int i = 0; i < state_num; i++){
  		alpha_init.push_back(pi[i] * b[i][seq_int[0][0]]); // FIX!!  LEFT 0 -> ONLY FIRST LINE
  	}
  	alpha.push_back(alpha_init);
  	
  	// Calculate values for each time period after initial alpha (forward algorithm)
  	for(int t = 1; t < time_period; t++){
  		// For each alpha[i][j]
  		for(int i = 0; i < state_num; i++){	
  			double tmp_sum = 0;
  			double tmp_b = 0;
  			for(int j = 0; j < state_num; j++){
  				tmp_sum += alpha[t - 1][j] * a[j][i];
  			}
  			alpha_value.push_back(tmp_sum * b[seq_int[0][t]][i]); // FIX!!	
  		}
  		// Store the row for time i to alpha
  		alpha.push_back(alpha_value);
  		alpha_value.clear();
  	}

  	// Create 2D vector beta - reverse //
  	vector< vector<double> > beta;
  	vector<double> beta_init, beta_value;
  	
  	// Set the initial values
  	for(int i = 0; i < state_num; i++){
  		beta_init.push_back(1);
  	}
  	beta.push_back(beta_init);

  	// Calculate values for each time period before initial beta (backward algorithm)
  	for(int t = 1; t < time_period; t++){
  		for(int i = 0; i < state_num; i++){
  			double tmp_sum = 0;
  			for(int j = 0; j < state_num; j++){
  				tmp_sum += a[i][j] * b[seq_int[0][time_period - t]][j] * beta[t-1][j];
  			}
  			beta_value.push_back(tmp_sum);
  		}
  		beta.push_back(beta_value);
  		beta_value.clear();
  	}
	reverse(beta.begin(),beta.end());

	// Create 2D vector gamma //
	vector< vector<double> > gamma;
	vector<double> gamma_value;

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
  				denom += alpha[t][i] * a[i][j] * b[seq_int[0][t + 1]][j] * beta[t + 1][j];
  			}	
  		}
  		for(int i = 0; i < state_num; i++){
  			for(int j = 0; j < state_num; j++){
  				epsilon_value.push_back(alpha[t][i] * a[i][j] * b[seq_int[0][t + 1]][j] * beta[t + 1][j] / denom);
  			}
  			epsilon_ij.push_back(epsilon_value);
  			epsilon_value.clear();
  		}
  		epsilon.push_back(epsilon_ij);
  		epsilon_ij.clear();
  	}

  	// gamma_one to store the sum of gamma ones for the pi update
	// gamma_t_min to store the sums of denominators for the a_ij update
	// gamma_t to store the sums of denominators for the b_jk update
	vector<double> gamma_one(state_num,0);
	vector<double> gamma_t_min(state_num,0);
	vector<double> gamma_t(state_num,0);

	// epsilon_t_min to store the sum of nominators for the a_ij update
	// gamma_t_each to store the sum of nominators for the b_jk update
	vector<vector<double> > epsilon_t_min(state_num,vector<double>(state_num,0));
	vector<vector<double> > gamma_t_each(observ_num,vector<double>(state_num,0));
	
  	// Accumulate gamma_one, epsilon_t_min, gamma_t_min, gamma_tk, gamma_t ...
  	for(int i = 0; i < state_num; i++){
  		gamma_one[i] = gamma[0][i];  // FIX?
  	}


  	for(int t = 0; t < time_period; t++){
  		// for epsilon_t_min and gamma_t_min
  		if(t != time_period - 1){
  			for(int i = 0; i < state_num; i++){
  				gamma_t_min[i] += gamma[t][i];
  				gamma_t[i] += gamma[t][i];
  				gamma_t_each[seq_int[0][t]][i] += gamma[t][i]; // FIX!!
  				for(int j = 0; j < state_num; j++){
  					epsilon_t_min[i][j] += epsilon[t][i][j];
  				}
  			}
  		}
  		else{
		// deal with last time period (gamma_t_each and gamma_t)
  			for(int i = 0; i < state_num; i++){
  				gamma_t[i] += gamma[t][i];
  				gamma_t_each[seq_int[0][t]][i] += gamma[t][i]; // FIX!!
  			}  			
  		}
  	}

  	// Update pi, a and b
  	for(int i = 0; i < state_num; i++){
  		
  		pi[i] = gamma_one[i];
  		for(int j = 0; j < state_num; j++){
  			a[i][j] = epsilon_t_min[i][j] / gamma_t_min[i];
  		}
  		for(int k = 0; k < observ_num; k++){
  			b[k][i] = gamma_t_each[k][i] / gamma_t[i];
  		}
  	}

  	

  	/* Test print 
  	// Printing pi, a and b
	printf("pi:\n");
  	for(int i = 0; i < state_num; i++){
  		printf("%f ", pi[i]);
  	}
  	printf("\n");

  	printf("a:\n");
  	for(int i = 0; i < state_num; i++){
  		for(int j = 0; j < state_num; j++){
  			printf("%f ", a[i][j]);	
  		}
  		printf("\n");
  	}
  	printf("\n");

  	printf("b:\n");
  	for(int i = 0; i < observ_num; i++){
  		for(int j = 0; j < state_num; j++){
  			printf("%f ", b[i][j]);	
  		}
  		printf("\n");
  	}
  	printf("\n");

  	/* Test print
  	// Printing accumulations
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

	for(int i = 0; i < observ_num; i++){
		printf("gamma %i:\n", i);
		for(int j = 0; j < state_num; j++){
			printf("%f ",gamma_t_each[i][j]);	
		}	
		printf("\n");		
	}
	/*
  	
  	
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



	
/*	
	HMM hmm_initial;
	loadHMM( &hmm_initial, "model_init.txt" );


	

	printf("%s\n", hmm_initial.model_name);
	printf("%i\n", hmm_initial.state_num);
	printf("%i\n", hmm_initial.observ_num);


	dumpHMM( stderr, &hmm_initial );
*/		

	return 0;
}
