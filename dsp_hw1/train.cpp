#include "hmm.h"
#include <math.h>
#include <iostream>
#include <fstream>
#include <string>
using namespace std;

int main()
{
/*
	HMM hmms[5];
	load_models( "modellist.txt", hmms, 5);
	dump_models( hmms, 5);
*/
	// Load seq_model
	vector<string> seq_model;
	string line;
  	ifstream myfile ("seq_model_01.txt");
  	if (myfile.is_open()){
    	while ( getline (myfile,line) ){
    		seq_model.push_back(line);
    	}
    	myfile.close();
  	}

  	



	
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
