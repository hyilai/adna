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

#define MIN_ADAPTER_LENGTH 7


class SmithWaterman;

int main (int argc, char** argv) {

	srand(time(NULL));

	if (argc != 4) {
		cerr << "Arguments: <result length> <minimum matching score> <output file>" << endl;
		exit(1); 
	}
	int length = atoi(argv[1]);
	int match_score = atoi(argv[2]);

	ofstream out(argv[3]);	// result file
	ifstream in1("testRead1.fastq");
	ifstream in2("testRead2.fastq");

	// containers
	vector<string> matched_seq;

	// analyze sw
	int j = 0;
	vector<double> total_elapsed_time;

	// clock start
	clock_t start = clock();

	int num_matched = 0;
	int total_derv = 0;
	int num_perfect = 0;

	string info1, sequence1, na1, quality1;
	string info2, sequence2, na2, quality2;
	while(getline(in1, info1) && getline(in2, info2)) {

		// read lines from in1
		getline(in1, sequence1);
		getline(in1, na1);
		getline(in1, quality1);

		// read lines from in2
		getline(in2, sequence2);
		getline(in2, na2);
		getline(in2, quality2);

		// cout << "Reading line " << j+1 << "..." << endl;

		string matched = "";
		int highest_score = 0;

		// start clock for sw algorithm
		clock_t start_time = clock();

		// check one direction
		SmithWaterman* sw = new SmithWaterman(sequence1, sequence2, quality1, quality2, match_score);
		int curr_high_score = sw->get_highest_score();
		if (highest_score < curr_high_score) {
			if (sw->match_reads()) {
				matched = sw->get_matched_string();
			}
		}
		delete sw;


		// end clock for sw algorithm
		clock_t end_time = clock();

		// for calculations
		total_elapsed_time.push_back(double(end_time - start_time) / CLOCKS_PER_SEC);

		// output file
		if (matched.length() > 0) {
			num_matched++;
		}
		total_derv += (abs((int)matched.length() - length));

		if (matched.length() == length) {
			num_perfect++;
		}

		out << matched << endl;
		
		//increment loop
		j++;
	}

	cout << "Total reads: " << j << endl;
	cout << "Number matched: " << num_matched << endl;
	cout << "Number of perfectly merged reads: " << num_perfect << endl;
	cout << "Average derivation from actual length: " << (double)total_derv/(double)j << endl;
	cout << "Average derivation from actual length for non-perfectly merged reads: " << (double)total_derv/(double)(j - num_perfect) << endl;

	// clock ends
	clock_t end = clock();
	double total_time_passed = double(end - start) / CLOCKS_PER_SEC;

	out.close();


	return 0;
}
