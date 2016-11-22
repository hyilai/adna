#include <cstdlib>
#include <cstdio>
#include <string.h>
#include <vector>
#include <iomanip>
#include <iostream>
#include <fstream>
#include <ctime>
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

	vector<int> matched_adapter;
	vector<string> matched_seq;
	vector<string> trimmed_extra;
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

	double total_elapsed_time = 0.0;
	int lines = 0;

	while(1) {
		if (file.peek() == EOF) {
			break;
		}

		cout << "Reading line " << j+1 << "..." << endl;
		j++;

		string line1, line2, line3, line4, trimmed, extra, matched;
		int match_adapt = -1;
		int highest_score = 0;

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

		clock_t start_time = clock();

		for (int i = 0; i < adapters.size(); i++) {

			SmithWaterman* sw = new SmithWaterman(line2,adapters[i], line4, "", 0);

			int curr_high_score = sw->get_highest_score();

			if (highest_score < curr_high_score) {
				highest_score = curr_high_score;
				trimmed = sw->trim_from_ending();
				extra = sw->get_trimmed();
				matched = sw->get_matched_string();
				match_adapt = i;
			}

			delete sw;
		}

		clock_t end_time = clock();

		// for calculations
		total_elapsed_time += double(end_time - start_time) / CLOCKS_PER_SEC;
		lines++;

		trimmed_seq.push_back(trimmed);
		trimmed_extra.push_back(extra);
		matched_seq.push_back(matched);
		matched_adapter.push_back(match_adapt);
	}

	file.close();


	string outfile("test_output.txt");
	ofstream out(outfile.c_str());
	ofstream out2("trimmed_extra.txt");

	for (int i = 0; i < trimmed_seq.size(); i++) {
		//out << trimmed_seq[i] << "\t\t" << trimmed_extra[i] << "\t\t" << matched_seq[i] << "\t\t" << adapters[matched_adapter[i] ] << "\n";
		//out << trimmed_seq[i] << "\t\t" << trimmed_extra[i] << "\t\t" << matched_seq[i] << "\n";
		out << trimmed_seq[i] << "\t\t" << trimmed_extra[i] << "\n";
		// out << trimmed_seq[i] << endl;
		// out2 << trimmed_extra[i] << "\t\t" matched_seq[i] << endl;
	} 
	
	out.close();
	out2.close();

	// for calculations
	double average_elapsed_time = total_elapsed_time / (double) lines;
	cout << "Average elapsed time for each line: " << average_elapsed_time << " s" << endl; 

	return 0;
}
