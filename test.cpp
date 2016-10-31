#include <cstdlib>
#include <cstdio>
#include <string.h>
#include <vector>
#include <iomanip>
#include <iostream>
#include <fstream>
#include "SmithWaterman.hpp"
// #include <pthread.h>

using namespace std;

class SmithWaterman;

int main (int argc, char** argv) {

	if (argc != 3) {
		cerr << "Arguments: <fastq_file> <adapter_file>" << endl;
		exit(1); 
	}

	ifstream file(argv[1]);
	ifstream adpt(argv[2]);

	vector<string> trimmed_seq;
	vector<string> adapters;

	// get the list of adapters
	string ad;
	while (1) {
		if (adpt.peek() == EOF) {
			break;
		}

		//discard first line
		getline(adpt, ad);

		//get adapter sequence in the second line
		getline(adpt, ad);
		adapters.push_back(ad);
	}

	int j = 0;

	while(1) {
		if (file.peek() == EOF) {
			break;
		}

		cout << "Reading line " << j++ << "..." << endl;

		string line1, line2, line3, line4, trimmed;

		// first line of read
		// information of the read
		getline(file, line1);

		// second line of read
		// sequence read
		getline(file, line2);

		// third line of read
		// skip this line
		getline(file, line3);


		// fourth line of read 
		// quality control
		getline(file, line4);

		for (int i = 0; i < adapters.size(); i++) {

			SmithWaterman* sw = new SmithWaterman(line2,adapters[i], line4, "", 0);

			// string temp = sw->trim_both_sides();
			string temp = sw->trim_from_ending();

			if (temp.length() < line2.length()) {
				trimmed = temp;
			}

			delete sw;
		}

		trimmed_seq.push_back(trimmed);
	}

	file.close();


	string outfile("test_output.txt");
	ofstream out(outfile.c_str());

	for (int i = 0; i < trimmed_seq.size(); i++) {
		out << trimmed_seq[i] << "\n";
	} 
	
	out.close();

	return 0;
}