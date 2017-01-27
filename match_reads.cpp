/************************************************************
 ****
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
#include <mpi.h>
#include "SmithWaterman.hpp"

using namespace std;

class SmithWaterman;

// containers


// void parse_file (int start, int end) {
void concat_reads () {



	//start clock
	clock_t start_time = clock();



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


	in.close();
	out.close();
	junk_out.close();

	return;
}


int main (int argc, char** argv) {

	if (argc < 5) {
		cerr << "Arguments: <fastq_file1> <fastq_file2> <output_file> <minimum match length>" << endl;
		exit(1); 
	}

	ifstream file1(argv[1]);
	ifstream file2(argv[2]);

	int minimum_read_length = atoi(argv[4]);
	int num_lines = -1;
	if (argc == 6) {
		num_lines = atoi(argv[5]);
	}


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

	concat_reads();


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



