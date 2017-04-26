#include <string.h>
#include <sstream>
#include <fstream>
#include <iostream>
#include <vector>
#include "utilities.hpp"
#include "global.hpp"

using namespace std;

// merge output files from all processes together
void merge_files (string outfile, int size) {

	// output name for the final output file
	stringstream ss;
	ss << outfile;
	string f = ss.str();
	ofstream out(f.c_str());


	for (int i = 0; i < size; i++) {

		// filename for thread i
		stringstream t;
		t << i << "_" << outfile;
		string curr_file = t.str();

		// merge
		ifstream in(curr_file.c_str());
		out << in.rdbuf();
		in.close();

		remove(curr_file.c_str());
	}

	out.close();
}


// merge output files from all workers together
void merge_worker_files (string outfile, int size) {

	// output name for the final output file
	stringstream ss;
	ss << outfile;
	string f = ss.str();
	ofstream out(f.c_str());


	for (int i = 1; i < size; i++) {

		// filename for thread i
		stringstream t;
		t << i << "_" << outfile;
		string curr_file = t.str();

		// merge
		ifstream in(curr_file.c_str());
		out << in.rdbuf();
		in.close();
		remove(curr_file.c_str());
	}

	out.close();
}


string get_file_name(string type, int file_num, string outfile) {
	stringstream fss;
	fss << type << "_" << file_num << "_" << outfile;
	string out_filename = fss.str();

	return out_filename.c_str();
}


string get_rank_file_name(string type, int file_num, int rank, string outfile) {
	stringstream fss;
	fss << rank << "_" << type << "_" << file_num << "_" << outfile;
	string out_filename = fss.str(); 

	return out_filename.c_str();
}


bool is_same_pair (string info1, string info2) {

	string info1_first, info1_second, info2_first, info2_second;

	size_t found = info1.find(" ");
	if (found) {
		info1_first = info1.substr(0, found);
		info1_second = info1.substr(found+1);
	}

	found = info2.find(" ");
	if (found) {
		info2_first = info2.substr(0, found);
		info2_second = info2.substr(found+1);
	}

	// cout << info1_first << endl;
	// cout << info2_first << endl;

	// check if the reads are of the same pair
	if (info1_first.compare(info2_first) == 0) {
		return true;
	} else {
		return false;
	}
}


bool read_from_fastq (ifstream &in, string &info, string &sequence, string &extra, string &quality) {

	// return false if the data is incomplete
	if (!getline(in, info)) return false; // get info
	if (!getline(in, sequence)) return false; // get read
	if (!getline(in, extra)) return false; // third line
	if (!getline(in, quality)) return false; // get the quality score

	return true;
}


void write_to_fastq (ofstream &file_stream, string info, string seq, string extra, string qual) {
	file_stream << info << endl
				<< seq << endl
				<< extra << endl
				<< qual << endl;
}


// print diagnostics after the program runs
void print_diagnostics (double elapsed_time, int num_reads, int discarded1, int discarded2, int num_concat, int num_final, int num_non_concat1, int num_non_concat2) {

	int elapsed_sec, elapsed_min, elapsed_hrs, elapsed_days;
	int seconds = (int) elapsed_time;
	int minutes = (seconds > 0) ? seconds / 60 : 0;
	int hours = (minutes > 0) ? minutes / 60 : 0;
	int days = (hours > 0) ? hours / 24 : 0;
	elapsed_sec = seconds % 60;
	elapsed_min = minutes % 60;
	elapsed_hrs = hours % 24;
	elapsed_days =  days;

	// print elapsed time
	cout << "Elapsed time: ";
	if (elapsed_days > 0) cout << elapsed_days << " days ";
	if (elapsed_hrs > 0) cout << elapsed_hrs << " hrs ";
	if (elapsed_min > 0) cout << elapsed_min << " min ";
	cout << elapsed_sec << " s" << endl;

	// print number of discarded reads
	cout << "Minimum read length after trimming: " << minimum_read_length << endl;
	cout << "Minimum matching length: " << minimum_match_length << endl;
	cout << "Total number of reads: " << num_reads << endl;
	cout << "Number of discarded reads in file 1 due to being trimmed too short: " << discarded1 << endl;
	cout << "Number of discarded reads in file 2 due to being trimmed too short: " << discarded2 << endl;
	cout << "Number of concatenated reads: " << num_concat << endl;
	cout << "NUmber of concatenated reads that passed quality check: " << num_final << endl;
	cout << "NUmber of non-concatenated reads from file 1 that passed quality check: " << num_non_concat1 << endl;
	cout << "NUmber of non-concatenated reads from file 2 that passed quality check: " << num_non_concat2 << endl;

}
