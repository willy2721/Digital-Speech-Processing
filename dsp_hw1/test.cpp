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
	if(argc != 4){
		printf("%s unable to execute due to wrong number of arguments\n", argv[0]);
	}
	
	char* model_list = argv[1];
	char* test_data = argv[2];
	char* out_result = argv[3];
	

	// Set up global variables//
	int time_period = 0;
	int state_num = 0;
	int observ_num = 0;
	int sample_num = 0;

	// Define a map to map chars to int //
	map<char,int> seq_map;
	seq_map['A'] = 0;
	seq_map['B'] = 1;
	seq_map['C'] = 2;
	seq_map['D'] = 3;
	seq_map['E'] = 4;
	seq_map['F'] = 5;

	// Load testing data //
	vector<string> test_seq;
	string line;
  	ifstream myfile(test_data);
  	if (myfile.is_open()){
    	while ( getline (myfile,line) ){
    		if(time_period == 0) time_period = line.length();
    		test_seq.push_back(line);
    		sample_num++;
    	}
    	myfile.close();
  	}

  	// Create integer representation of test_seq (A = 0 ... F = 5) // 
  	vector< vector<int> > test_int;
  	for(int i = 0; i < test_seq.size(); i++){
  		vector<int> tmp;
  		for(int j = 0; j < test_seq[i].size(); j++){
  			tmp.push_back(seq_map.find(test_seq[i][j])->second);
  		}
  		test_int.push_back(tmp);
  	}

  	HMM hmms[5];
	load_models(model_list, hmms, 5);
	state_num = hmms[0].state_num;
	observ_num = hmms[0].observ_num;

	//dump_models( hmms, 5);
	// vector to store answer for every sample
	vector<int> ans;
	vector< vector<double> > delta;
	vector<double > delta_value;

	for(int n = 0; n < sample_num; n++){
		double max_prob = 0;
		int max_model = -1;
		
		// Try the first hmm for the first sample
		for(int i = 0; i < 5; i++){
			// Record the largest probability for each model	
			double max_delta = 0;
			for(int t = 0; t < time_period; t++){
				for(int j = 0; j < state_num; j++){
					double max_delta_a = 0;
					if(t == 0){
						delta_value.push_back(hmms[i].initial[j] * hmms[i].observation[test_int[n][t]][j]);	
					}
					
					//printf("%f",hmms[0].initial[j] * hmms[0].observation[test_int[0][0]][j]);
					else {
						for(int k = 0; k < state_num; k++){
							if(delta[t - 1][k] * hmms[i].transition[k][j] > max_delta_a){
								max_delta_a = delta[t - 1][k] * hmms[i].transition[k][j];
							}
						}
						delta_value.push_back(max_delta_a * hmms[i].observation[test_int[n][t]][j]);
					}	
				}
				delta.push_back(delta_value);
				delta_value.clear();
				if(t == time_period - 1){
					for(int s = 0; s < state_num; s++){
						if(delta[t][s] > max_delta)
							max_delta = delta[t][s];
					}
				}
			}
			
			if(max_delta > max_prob){
				max_prob = max_delta;
				max_model = i;
			}

			delta.clear();
			
		}
		ans.push_back(max_model);
	}

	ofstream outfile;
	outfile.open(out_result);
	string output = "";
	for(int i = 0; i < sample_num; i++){
		output += "model_0" + patch::to_string(ans[i] + 1) + ".txt\n";
	}
	
	outfile << output;
	outfile.close();

	//for(int i = 0; i < ans.size(); i++){
	//	printf("%i ",ans[i]);	
	//}
	

	/*
	for(int i = 0; i < time_period; i++){
s		for(int j = 0; j < state_num; j++){
			printf("%f ",delta[i][j]);
		}
		printf("\n");
	}
	*/







	return 0;
}
