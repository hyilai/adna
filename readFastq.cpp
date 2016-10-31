/************************************************************
 ****
 ****	This portion reads in a fastq file and fa file
 ****	and trims the adapter sequences from the reads
 ****	in the fastq file. It outputs an output file
 ****	with the trimmed reads along with their quality
 ****	control strings.
 ****
 ************************************************************/



#include <cstdlib>
#include <cstdio>
#include <string.h>
#include <vector>
#include <iomanip>
#include <iostream>
#include <fstream>
#include "SmithWaterman.hpp"

using namespace std;

class SmithWaterman;

int main (int argc, char** argv) {

	if (argc != 4) {
		cerr << "Arguments: <fastq_file> <adapter_file> <output_file>" << endl;
		exit(1); 
	}

	ifstream file(argv[1]);
	ifstream adpt(argv[2]);

	vector< vector<string> > trimmed_seq;
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


	ofstream out(argv[3]);

	int j = 0;

	// read fastq file
	while(1) {
		if (file.peek() == EOF) {
			break;
		}

		cout << "Reading line " << j++ << "..." << endl;
		string line1, line2, line3, line4, trimmed, quality;

		// first line of read
		// information of the read
		getline(file, line1);

		// second line of read
		// sequence read
		getline(file, line2);

		// third line of read
		getline(file, line3);

		// fourth line of read 
		// quality control
		getline(file, line4);

		for (int i = 0; i < adapters.size(); i++) {

			SmithWaterman* sw = new SmithWaterman(line2,adapters[i], line4, "", 0);

			// string temp = sw->trim_both_sides();	// trim read

			if (temp.length() > 0 && temp.length() < line2.length()) {
				trimmed = temp;
				quality = sw->get_quality1();	// get trimmed quality control string
			}

			trimmed_seq[j].push_back(temp);

			delete sw;
		}

		// output new data
		out << line1 << endl
			<< trimmed << endl
			<< line3 << endl
			<< quality << endl;
	}

	file.close();
	out.close();

	return 0;
}