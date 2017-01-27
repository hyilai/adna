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
#include <ctime>
// #include <mpi.h>
#include "SmithWaterman.hpp"

using namespace std;

// universal sequence for tru seq adapters
const string TRU_SEQ_ADAPTER = "GATCGGAAGAGCACACGTCTGAACTCCAGTCAC";

class SmithWaterman;

// containers
vector<vector<string> > file_line;
vector< vector<string> > trimmed_seq;
vector<string> adapters;

// void parse_file (int start, int end) {
void trim_reads (char* infile, char* outfile, int min_length) {

	int discarded = 0;

	ifstream in(infile);
	ofstream out(outfile);
	string junk_filename = "junk_" + string(outfile);
	ofstream junk_out(junk_filename.c_str());

	// trim fastq file
	// read 4 lines in a roll (4 lines = 1 set of data)
	// for (int i = start; i < end; i++) {

	//start clock
	clock_t start_time = clock();

	int j = 0;
	string info;
	while (getline(in, info)) {

		if (j % 100 == 0) {
			cout << "Reading line " << j+1 << "..." << endl;
		}
		
		j++;
		
		string sequence, extra, quality;
		string read, trimmed_junk, junk_quality, new_quality;
		int highest_score = 0;

		// read lines from the input file
		// getline(in, info); // get fastq read info
		getline(in, sequence); // get read info
		getline(in, extra); // N/A
		getline(in, quality); // get the quality scores
		

		// for (int i = 0; i < adapters.size(); i++) {

			SmithWaterman* sw = new SmithWaterman(sequence, TRU_SEQ_ADAPTER, quality, "", 0);

			int curr_high_score = sw->get_highest_score();

			if (highest_score < curr_high_score) {
				highest_score = curr_high_score;
				read = sw->trim_from_ending();
				new_quality = sw->get_quality1();
				trimmed_junk = sw->get_trimmed();
				junk_quality = sw->get_trimmed_quality();
			}

			delete sw;
		// }

		if (read.length() >= min_length) {
			// write to output file
			out << info << endl
				<< read << endl
				<< extra << endl
				<< new_quality << endl;

			// write to junk output file
			junk_out << info << endl
					 << trimmed_junk << endl
					 << extra << endl
					 << junk_quality << endl;
		} else {
			discarded++;
		}

	}

	// end timer
	clock_t end_time = clock();

	// find elapsed time
	double elapsed_time = (end_time - start_time) / (double) CLOCKS_PER_SEC;
	int elapsed_sec, elapsed_min, elapsed_hrs, elapsed_days;
	elapsed_sec = (int) elapsed_time;
	elapsed_min = (elapsed_sec > 0) ? elapsed_sec / 60 : 0;
	elapsed_hrs = (elapsed_min > 0) ? elapsed_min / 60 : 0;
	elapsed_days = (elapsed_hrs > 0) ? elapsed_hrs / 24 : 0;

	// print elapsed time
	cout << "Elapsed time: ";
	if (elapsed_days > 0) cout << elapsed_days << " days ";
	if (elapsed_hrs > 0) cout << elapsed_hrs << " hrs ";
	if (elapsed_min > 0) cout << elapsed_min << " min ";
	cout << elapsed_sec << " s" << endl;

	// print number of discarded reads
	cout << "Total number of reads: " << j << endl;
	cout << "Number of discarded reads due to being trimmed too short: " << discarded << endl;

	in.close();
	out.close();
	junk_out.close();

	return;
}


int main (int argc, char** argv) {

	if (argc < 5) {
		cerr << "Arguments: <fastq_file> <adapter_file> <output_file> <minimum length of trimmed read> <number of lines (optional)>" << endl;
		exit(1); 
	}

	// ifstream file(argv[1]);
	ifstream adpt(argv[2]);

	int minimum_read_length = atoi(argv[4]);
	int num_lines = -1;
	if (argc == 6) {
		num_lines = atoi(argv[5]);
	}

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


	//Match 2 reads and put them together
	// string info;
	// while (getline(file, info)) {

	// 	string sequence, line3, quality, matched_string;
	// 	getline(file, sequence);
	// 	getline(file, line3);
	// 	getline(file, quality);

	// 	for (int i = 0; i < adapters.size(); i++) {
	// 		SmithWaterman* sw = new SmithWaterman(sequence, adapters[i], quality, "", minimum_match_score);


	// 		bool is_Matched = false;
	// 		int curr_high_score = sw->get_highest_score();

	// 		if (highest_score < curr_high_score) {
	// 			highest_score = curr_high_score;
	// 			is_Matched = sw->match_reads();
	// 			if(is_Matched) {
	// 				matched_string = sw->get_matched_string();
	// 			}
	// 			else {
	// 				matched_string = "";
	// 			}
	// 		}

	// 		delete sw;
	// 	}
	// }
	

	// get all file lines
	// string line;
	// while (getline(file,line)) {
	// 	vector<string> line_group;
	// 	line_group.push_back(line);	// first line; information of the read
	// 	getline(file, line);
	// 	line_group.push_back(line);	// second line; sequence read
	// 	getline(file, line);
	// 	line_group.push_back(line);	// third line
	// 	getline(file, line);
	// 	line_group.push_back(line);	// fourth line; quality control
	// 	file_lines.push_back(line_group);
	// }


	// ofstream out(argv[3]);
	// string junk_filename = "junk" + string(argv[3]);
	// ofstream out(junk_filename.c_str());

	trim_reads(argv[1], argv[3], minimum_read_length);


	// set up MPI stuff
	// int comm_sz;
	// int my_rank;

	// MPI_Init(NULL,NULL);
	// MPI_comm_size(MPI_COMM_WORLD, &comm_sz);
	// MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

	// if (my_rank != 0) {
	// 	//MPI_Send();
	// } else {
	// 	//MPI_Recv();
	// }

	// MPI_Finalize();

	return 0;
}



