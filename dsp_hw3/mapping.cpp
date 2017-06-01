#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <map>
#include <algorithm>
#include <boost/algorithm/string/split.hpp>                                      
#include <boost/algorithm/string.hpp>
using namespace std;


int main(int argc, char *argv[])
{
	char* in_map = argv[1];
	char* out_map = argv[2];
	map<string,vector<string> > mymap;
	string line, big_c, zhu_yin, zhu;
	vector<string> results;												// Keep results of zhu_yin split
	ifstream myfile (in_map);
  	if (myfile.is_open()){
    	while ( getline (myfile,line) ){
    		big_c = line.substr(0,2);									// Extract chinese character
    		zhu_yin = line.substr(3);									// Extract zhu_yin
			boost::split(results, zhu_yin, boost::is_any_of("/"));		// C++ equivalent of python string split
			for(int i = 0; i < results.size(); i++){
				zhu = results[i].substr(0,2);							// Extract first two characters
				// Add only if unique
				if(find(mymap[zhu].begin(), mymap[zhu].end(), big_c) == mymap[zhu].end())
					mymap[zhu].push_back(big_c);
			}
			mymap[big_c].push_back(big_c);								// Also add character based mapping
    	}
    	myfile.close();
  	}

  	map<string, vector<string> >::iterator iter;
  	ofstream outfile;
	outfile.open(out_map);
	string output = "";
	for (iter = mymap.begin(); iter != mymap.end(); iter++) {
  		output += iter->first;
  		output += "      ";
  		for(int i = 0; i < iter->second.size(); i++){
  			output += iter->second[i];
  			if(i < iter->second.size() - 1) output += " ";
  		}
  		output += "\n";
  	}
	outfile << output << endl;
	outfile.close();

	return 0;
}
