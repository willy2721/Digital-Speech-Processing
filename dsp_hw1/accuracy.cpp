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



int main(int argc, char *argv[])
{
	// Reads two text files and compares the accuracy
	if(argc != 3){
		printf("%s unable to execute due to wrong number of arguments\n", argv[0]);
	}
	
	char* result = argv[1];
	char* answer = argv[2];
	vector<int > res;
	vector<int > ans;
	double ac = 0;
	double num = 0;
	double accuracy = 0;

	string line;
  	ifstream myresult(result);
  	if (myresult.is_open()){
    	while ( getline (myresult,line) ){
    		num++;
    		if(line == "model_01.txt")
    			res.push_back(1);
    		else if(line == "model_02.txt")
    			res.push_back(2);
    		else if(line == "model_03.txt")
    			res.push_back(3);
    		else if(line == "model_04.txt")
    			res.push_back(4);
    		else if(line == "model_05.txt")
    			res.push_back(5);

    	}
    	myresult.close();
  	}

  	ifstream myanswer(answer);
  	if (myanswer.is_open()){
    	while ( getline (myanswer,line) ){
    		if(line == "model_01.txt")
    			ans.push_back(1);
    		else if(line == "model_02.txt")
    			ans.push_back(2);
    		else if(line == "model_03.txt")
    			ans.push_back(3);
    		else if(line == "model_04.txt")
    			ans.push_back(4);
    		else if(line == "model_05.txt")
    			ans.push_back(5);

    	}
    	myanswer.close();
  	}

  	for(int i = 0; i < ans.size(); i++){
  		if(res[i] == ans[i])
  			ac++;
  	}
  	accuracy = ac / num;
  	printf("Correct answers = %f\n Total answers = %f\n Accuracy = %f\n",ac,num,accuracy);
	//char* out_result = argv[3];


	return 0;
}