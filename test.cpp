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

char get_random_base() {
	int random = rand() % 20 + 1;
	if(random < 6) {
		return 'A';
	} 
	else if(random < 11) {
		return 'G';
	} 
	else if(random < 16) {
		return 'T';
	}
	else {
		return 'C';
	}
}


//Generates DNA of random length between 0 and range
string make_extra_DNA(int range) {
	string DNA;

	int i;
	stringstream ss;
	for(i=0;i<range;i++) {
		ss << get_random_base();
	}

	ss >> DNA;

	ss.str("");
	ss.clear();

	return DNA;
}


string make_read (int read_length, string adapter, int &adapter_length, string &adapter_fragment, string &extra_fragment) {

	// generate random read sequence
	string sequence = make_extra_DNA(read_length);

	// cut a chunk of the adapter out
	// length of adapter fragment is between 7 to length of adapter
	// randomly add bases to the end of adapter
	adapter_length = rand() % (adapter.length() - MIN_ADAPTER_LENGTH) + MIN_ADAPTER_LENGTH;
	adapter_fragment = adapter.substr(0,adapter_length);

	//Generates extra DNA for the ends of the adapters
	//The length of the extra DNA will be between 0 - 3/4 of adapter fragment
	int num_extra = rand() % ((adapter_length*3)/4);
	if (num_extra) {
		extra_fragment = make_extra_DNA(num_extra);
	} else {
		extra_fragment = "";
	}

	sequence += adapter_fragment + extra_fragment;

	return sequence;
}

class SmithWaterman;

int main (int argc, char** argv) {

	srand(time(NULL));

	if (argc != 6) {
		cerr << "Arguments: <number of fragments> <read length> <minimum length of cut read> <adapter file> <output file>" << endl;
		exit(1); 
	}

	int num_lines = atoi(argv[1]);
	int read_length = atoi(argv[2]);
	int min_cut_length = atoi(argv[3]);

	ifstream adpt(argv[4]);
	if (!adpt) {
		cerr << "Error: cannot open adapter file" << endl;
		exit(1);
	}

	ofstream out(argv[5]);	// result file
	ofstream testfile("test_sequences.txt");	// file containing test sequences
	ofstream stat_file("statistics.csv");		// file dump for a bunch of numbers

	// containers
	vector<int> matched_adapter;
	vector<int> actual_adapter;
	vector<string> trimmed_extra;
	vector<string> trimmed_seq;
	vector<string> matched_seq;
	vector<string> actual_seq;
	vector<string> actual_extra;
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
	adpt.close();

	// analyze sw
	int j = 0;
	vector<double> total_elapsed_time;
	stat_file << ",time_passed,actual_read_length,trimmed_seq_length,actual_adapter_frag_length,actual_extra_length,discarded_seq_length,matched_seq_length,actual_adapter_num,matched_adapter_num" << endl;

	// clock start
	clock_t start = clock();

	while(j < num_lines) {
		cout << "Reading line " << j+1 << "..." << endl;
		
		/* make sequence */
		int ad_num = rand() % adapters.size();
		string adapter = adapters[ad_num];
		int ad_frag_length;
		string ad_frag, extra;
		string sequence = make_read(read_length, adapter, ad_frag_length, ad_frag, extra);
		string quality = make_extra_DNA(sequence.length());

		// add actual sequence structure to test file
		testfile << sequence.substr(0,read_length) << " " << ad_frag << " " << extra << endl;
		actual_seq.push_back(sequence.substr(0, read_length));
		actual_extra.push_back(ad_frag + extra);
		actual_adapter.push_back(ad_num);

		string trimmed, new_extra, matched;
		int matched_adapter_num;
		int highest_score = 0;

		// start clock for sw algorithm
		clock_t start_time = clock();

		for (int i = 0; i < adapters.size(); i++) {

			SmithWaterman* sw = new SmithWaterman(sequence, adapters[i], quality, "", 0);

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

		// for calculations
		total_elapsed_time.push_back(double(end_time - start_time) / CLOCKS_PER_SEC);

		trimmed_seq.push_back(trimmed);
		trimmed_extra.push_back(new_extra);
		matched_seq.push_back(matched);
		matched_adapter.push_back(matched_adapter_num);

		// output file
		out << trimmed << " " << new_extra << endl;

		// add to stat file
		//stat_file << ",time_passed,actual_read_length,trimmed_seq_length,actual_adapter_frag_length,actual_extra_length,discarded_seq_length,matched_seq_length,actual_adapter_num,matched_adapter_num" << endl;
		stat_file << j << ","
				  << total_elapsed_time[j] << "," 		// time_passed
				  << read_length << ","					// actual_read_length
				  << trimmed.length() << ","			// trimmed_seq_length
				  << ad_frag.length() << ","			// actual_adapter_frag_length
				  << extra.length() << ","				// actual_extra_length
				  << new_extra.length() << ","			// discared_seq_length
				  << matched.length() << ","			// matched_seq_length
				  << ad_num << ","						// actual_adapter_num
				  << matched_adapter_num << endl;		// matched_adapte_num

		//increment loop
		j++;
	}

	// clock ends
	clock_t end = clock();
	double total_time_passed = double(end - start) / CLOCKS_PER_SEC;

	testfile.close();
	out.close();


	// for calculations
	// ofstream result_file("result.txt");
	cout << endl;
	cout << "Result: " << endl;
	// out << "Result: " << endl;
	double total_calc_time = 0;
	int num_correct_adpt_match = 0;
	int num_correct_trimmed = 0;
	int num_total_trimmed_length = 0;
	int num_total_read_length = 0;
	vector<double> actual_trimmed_difference;
	double biggest_diff = 0.0;
	int biggest_diff_line = 0;
	for (int i = 0; i < num_lines; i++) {
		total_calc_time += total_elapsed_time[i];
		if (actual_adapter[i] == matched_adapter[i]) {
			num_correct_adpt_match++;
		}
		if (read_length == trimmed_seq[i].length() ) {
			num_correct_trimmed++;
		}
		num_total_trimmed_length += trimmed_seq[i].length();
		num_total_read_length += read_length;
		actual_trimmed_difference.push_back(trimmed_seq[i].length() / (double) read_length);
		if (abs(1.00 - actual_trimmed_difference[i]) > biggest_diff) {
			biggest_diff = actual_trimmed_difference[i];
			biggest_diff_line = i+1;
		}
	}
	double average_elapsed_time = total_calc_time / (double) num_lines;
	double average_trimmed_length = num_total_trimmed_length / (double) num_lines;
	cout << "Number of sequences: " << num_lines << endl;
	cout << "Actual read length: " << read_length << endl;
	cout << "Total time: " << total_calc_time << " s" << endl;
	cout << "Average Elapsed Time: " << average_elapsed_time << " s" << endl;
	// out << "Average Elapsed Time: " << average_elapsed_time << endl;
	cout << "Number of correctly matched adapters: " << num_correct_adpt_match << " out of " << num_lines << endl;
	cout << "Number of correctly trimmed sequences: " << num_correct_trimmed << " out of " << num_lines << endl;
	cout << "Average length of trimmed reads: " << average_trimmed_length << " of " << read_length << endl;
	cout << endl;
	if (biggest_diff_line != 0) {
		cout << "Sequence with worst trimmed result: line " << biggest_diff_line << endl;
		cout << "  Length of sequence after trimming: \t" << trimmed_seq[biggest_diff_line-1].length() << " of " << read_length << endl;
		cout << "  Actual read sequence: \t" << actual_seq[biggest_diff_line-1] << endl;
		cout << "  Sequence after trimming: \t" << trimmed_seq[biggest_diff_line-1] << endl;
		cout << "  Actual adapter: \t" << adapters[actual_adapter[biggest_diff_line-1] ] << endl;
		cout << "  Matched adapter: \t" << adapters[matched_adapter[biggest_diff_line-1] ] << endl;
		cout << "  Actual adapter fragment + extra sequence: \t" << " " << actual_extra[biggest_diff_line-1] << endl;
		cout << "  Matched sequence within original string: \t" << matched_seq[biggest_diff_line-1] << endl;
	} else {
		cout << "All lines are trimmed correctly" << endl;
	}
	

	return 0;
}
