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
	vector<double> pi;
	vector<vector<double> > a, b;
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
  	ifstream myfile;
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
	load_models( "modellist.txt", hmms, 5);
	dump_models( hmms, 5);
	
	// vector to store answer for every sample
	vector<int> ans;
	
	vector< vector<double> > delta;
	for(int n = 0; n < sample_num; n++){
		for(int i = 0; i < sizeof(hmms); i++){
			double max_prob = 0;
			int max_model;
			
		}
	}


	return 0;
}
