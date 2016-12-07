#include <cstdlib>
#include <cstdio>
#include <string.h>
#include <vector>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <fstream>
#include <ctime>
#include "SmithWaterman.hpp"
// #include <pthread.h>

using namespace std;


class SmithWaterman;

int main (int argc, char** argv) {

	srand(time(NULL));

	if (argc != 5) {
		cerr << "Arguments: <number of lines> <input file> <adapter file> <output file>" << endl;
		exit(1); 
	}

	int num_lines = atoi(argv[1]);

	ifstream in(argv[2]); // Input file
	ifstream in_adapter(argv[3]); // adapter file
	ofstream out(argv[4]);	// result file


	// containers
	vector<string> trimmed_extra;
	vector<string> trimmed_seq;
	vector<string> matched_seq;
	vector<string> actual_seq;
	vector<string> actual_extra;
	vector<string> adapters;



	string string_adapter;
	
	while (getline(in_adapter, string_adapter) ) {
		getline(in_adapter, string_adapter);
		adapters.push_back(string_adapter);
	}
	

	// clock start
	clock_t start = clock();

	int j = 0;
	while(j < num_lines) {
		cout << "Reading line " << j+1 << "..." << endl;
		
		string trimmed, new_extra, matched, info, read, umm, quality;
		int matched_adapter_num;
		int highest_score = 0;



		getline(in, info); // get fastq read info
		getline(in, read); // get read info
		getline(in, umm); // nothing 
		getline(in, quality); // get the quality scores
		

		// start clock for sw algorithm
		clock_t start_time = clock();

		for (int i = 0; i < adapters.size(); i++) {

			SmithWaterman* sw = new SmithWaterman(read, adapters[i], quality, "", 0);

			int curr_high_score = sw->get_highest_score();

			if (highest_score < curr_high_score) {
				highest_score = curr_high_score;
				trimmed = sw->trim_from_ending();
				new_extra = sw->get_trimmed();
				matched = sw->get_matched_string();
				matched_adapter_num = i;
			}

			delete sw;
		}


		// end clock for sw algorithm
		clock_t end_time = clock();


		// output file
		out << info << endl;
		out << trimmed << endl;
		out << umm << endl;
		out << quality << endl;

		//increment loop
		j++;
	}

	// clock ends
	clock_t end = clock();
	double total_time_passed = double(end - start) / CLOCKS_PER_SEC;

	
	out.close();


	
	return 0;
}
