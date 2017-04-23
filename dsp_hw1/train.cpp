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

	// Load initial HMM //
	HMM hmm_initial;
	loadHMM( &hmm_initial, "model_init.txt" );

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

  	// Create 2D vector alpha //
  	vector< vector<double> > alpha;
  	vector<double> alpha_init, alpha_value;

  	// Calculate the initial values
  	for(int i = 0; i < hmm_initial.state_num; i++){
  		double tmp_pi = hmm_initial.initial[i];
  		double tmp_b = hmm_initial.observation[i][seq_int[0][0]];
  		alpha_init.push_back(tmp_pi * tmp_b); // FIX!!  LEFT 0 -> ONLY FIRST LINE
  	}
  	alpha.push_back(alpha_init);
  	
  	// Calculate values for each time period after initial alpha (forward algorithm)
  	for(int i = 1; i < time_period; i++){
  		// For each alpha[i][j]
  		for(int j = 0; j < hmm_initial.state_num; j++){	
  			double tmp_sum = 0;
  			double tmp_b = 0;
  			for(int k = 0; k < hmm_initial.state_num; k++){
  				tmp_sum += alpha[i - 1][k] * hmm_initial.transition[k][j];
  			}
  			tmp_b = hmm_initial.observation[seq_int[0][i-1]][j]; // FIX
  			alpha_value.push_back(tmp_sum * tmp_b);	
  		}
  		// Store the row for time i to alpha
  		alpha.push_back(alpha_value);
  		alpha_value.clear();
  	}

  	// Create 2D vector beta - reverse //
  	vector< vector<double> > beta;
  	vector<double> beta_init, beta_value;
  	
  	// Set the initial values
  	for(int i = 0; i < hmm_initial.state_num; i++){
  		beta_init.push_back(1);
  	}
  	beta.push_back(beta_init);

  	// Calculate values for each time period before initial beta (backward algorithm)
  	for(int i = 1; i < time_period; i++){
  		for(int j = 0; j < hmm_initial.state_num; j++){
  			double tmp_sum = 0;
  			for(int k = 0; k < hmm_initial.state_num; k++){
  				tmp_sum += hmm_initial.transition[j][k] * hmm_initial.observation[seq_int[0][time_period - i]][k] * beta[i-1][k];
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
  	for(int i = 0; i < time_period; i++){
  		double denom = 0;
  		for(int j = 0; j < hmm_initial.state_num; j++){
  			denom += alpha[i][j] * beta[i][j];
  		}
  		for(int j = 0; j < hmm_initial.state_num; j++){
  			gamma_value.push_back(alpha[i][j] * beta[i][j] / denom);
  		}
  		gamma.push_back(gamma_value);
  		gamma_value.clear();
  	}

  	for(int i = 0; i < 5; i++){
  		for(int j = 0; j < 6; j++){
  			printf("%f ",gamma[i][j]);
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
  	for(int i = 0; i < 5; i++){
  		for(int j = 0; j < 6; j++){
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
