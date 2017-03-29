#include <string.h>
#include <string>
#include <fstream>
#include <vector>
#include "global.hpp"

using namespace std;


vector<string> adapters;

void set_adapters (char* c_filename) {
	// get the list of adapters
	ifstream adpt(c_filename);
	string ad;
	vector<string> adapters;
	while (1) {
		// if (adpt.peek() == EOF) {
		// 	break;
		// }

		//discard first line
		if (!getline(adpt, ad)) break;

		//get adapter sequence in the second line
		getline(adpt, ad);
		adapters.push_back(ad);
	}
	adpt.close();
}