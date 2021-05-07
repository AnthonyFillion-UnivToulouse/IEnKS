#include <iostream>
#include <string>
#include <map>
#include <vector>

#include "armadillo"
#include "cnpy.h"
#include "algos.h"

using namespace arma;
using namespace std;

int main (int argc, char *argv[]) {
	ifstream file(argv[1]);
	stringstream streamline;
	string line,key,value;
	map<string, string> icfg;
	while (getline(file, line))
	{
		streamline.clear();
		streamline.str(line);
		getline(streamline,key,'=');
		key.erase(remove(key.begin(), key.end(), ' '), key.end());
		if (streamline.rdbuf()->in_avail() != 0){
			streamline >> value;
			value.erase(remove(value.begin(), value.end(), ' '), value.end());
			icfg[key] = value;
		}
	}

	arma_rng::set_seed_random();

	//Algo
	map<string,double> ocfg;
	ocfg = algo_selector(icfg);
	string npy_file;
	const unsigned int shape[] = {1};
	for(map<string,double>::iterator it = ocfg.begin(); it != ocfg.end(); ++it) {
		npy_file = (icfg["directory"] +
				string("/data/") +
				icfg["name"] +
				string("_return=") +
				it->first +
				string(".npy")
				).c_str();
		remove(npy_file.c_str());
		cnpy::npy_save(npy_file,&(it->second),shape,1,"w");
	}
}
