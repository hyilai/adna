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
#include <string>
#include <vector>
#include <iomanip>
#include <iostream>
#include <fstream>
#include <sstream>
#include <ctime>
#include <unordered_map>
#include <mpi.h>
// #include <boost/iostreams/filtering_stream.hpp>
// #include <boost/iostreams/copy.hpp>
// #include <boost/iostreams/filter/gzip.hpp>
#include "memory_usage.hpp"
#include "steps.hpp"
#include "utilities.hpp"
#include "global.hpp"
#include "hash_table.hpp"
#include "MPI_readFastq.hpp"

using namespace std;

// tags for MPI messages
const int NUM_DISCARDED1_TAG = 0;
const int NUM_DISCARDED2_TAG = 1;
const int NUM_CONCAT_TAG = 2;
const int NUM_FINAL_TAG = 3;
const int NUM_NON_CONCAT1_TAG = 4;
const int NUM_NON_CONCAT2_TAG = 5;


// numbers for diagnostics
int num_reads1, num_reads2, discarded1, discarded2, num_concat, num_final, num_non_concat1, num_non_concat2;


int gather_number (int original_value, int source, int tag) {
	MPI_Status status;
	int temp;
	MPI_Recv(&temp, 1, MPI_INT, source, tag, MPI_COMM_WORLD, &status);
	return original_value + temp;
}


// single threaded read processing function
void process_reads (char* infile1, char* infile2, char* outfile, int rank, int size, vector<string> adapters, bool debug) {

	num_reads1 = 0; num_reads2 = 0; discarded1 = 0; discarded2 = 0; num_concat = 0; num_final = 0; num_non_concat1 = 0; num_non_concat2 = 0;

	//start clock
	clock_t start_time = clock();


	// diagnostic files
	string out_filename1 = get_rank_file_name(string("trimmed"), 1, rank, string(outfile));
	ofstream out1(out_filename1);
	string junk_filename1 = get_rank_file_name(string("junk"), 1, rank, string(outfile));
	ofstream junk_out1(junk_filename1);
	string discarded_filename1 = get_rank_file_name(string("discarded"), 1, rank, string(outfile));
	ofstream discarded_out1(discarded_filename1);
	string junk_discarded_filename1 = get_rank_file_name(string("junk_discarded"), 1, rank, string(outfile));
	ofstream junk_discarded_out1(junk_discarded_filename1);


	string out_filename2 = get_rank_file_name(string("trimmed"), 2, rank, string(outfile));
	ofstream out2(out_filename2);
	string junk_filename2 = get_rank_file_name(string("junk"), 2, rank, string(outfile));
	ofstream junk_out2(junk_filename2);
	string discarded_filename2 = get_rank_file_name(string("discarded"), 2, rank, string(outfile));
	ofstream discarded_out2(discarded_filename2);
	string junk_discarded_filename2 = get_rank_file_name(string("junk_discarded"), 2, rank, string(outfile));
	ofstream junk_discarded_out2(junk_discarded_filename2);


	// file for non-concatenated reads
	string non_concat_filename1 = get_rank_file_name(string("non_concat"), 1, rank, string(outfile));
	ofstream non_concat_out1(non_concat_filename1.c_str());
	string non_concat_filename2 = get_rank_file_name(string("non_concat"), 2, rank, string(outfile));
	ofstream non_concat_out2(non_concat_filename2.c_str());

	// file for concatenated reads
	stringstream ss;
	ss << rank << "_concat_" << string(outfile);
	string concat_filename = ss.str();
	ofstream concat_out(concat_filename.c_str());


	ifstream in1(infile1);
	ifstream in2(infile2);
	// open and read from gzipped file
	// ifstream file2(infile2, ios_base::in | ios_base::binary);
	// boost::iostreams::filtering_istream in2;
	// in2.push(boost::iostreams::gzip_decompressor());
	// in2.push(file2);

	bool has_more1, has_more2;


	// string info1, info2;
	// while (getline(in1, info1) && getline(in2, info2)) {
	while(1) {
		
		string info1, sequence1, extra1, quality1;
		string read1, trimmed_junk1, junk_quality1, new_quality1;
		string info2, sequence2, extra2, quality2;
		string read2, trimmed_junk2, junk_quality2, new_quality2;

		has_more1 = read_from_fastq(in1, info1, sequence1, extra1, quality1);
		has_more2 = read_from_fastq(in2, info2, sequence2, extra2, quality2);

		if (!has_more1 || !has_more2) break;

		// don't process this set of data if it's not for this process rank
		num_reads1++;
		num_reads2++;
		if ((num_reads1-1) % size != rank) continue;


		if (num_reads1 % LINEBLOCKS == 0) {
			cout << "Parsing read " << num_reads1 << "..." << endl;
		}

		if (!is_same_pair(info1, info2)) {
			cout << "\tWarning: read " << num_reads1 << " from file1 and file2 are of not the same pair; skipping" << endl;
			continue;
		}

		//strip T's
		strip_t(sequence1, quality1);
		strip_t(sequence2, quality2);
		
		// trim
		read1 = trim_read (sequence1, quality1, new_quality1, trimmed_junk1, junk_quality1, adapters);
		read2 = trim_read (sequence2, quality2, new_quality2, trimmed_junk2, junk_quality2, adapters);

		// flags for trimming
		bool read1_good,read2_good;

		/* write to debugging files */
		// read1
		if (read1.length() >= minimum_read_length) {
			// write to output file
			if (debug) write_to_fastq(out1, info1, read1, extra1, new_quality1);

			// write to junk output file
			if (debug) write_to_fastq(junk_out1, info1, trimmed_junk1, extra1, junk_quality1);

			read1_good = true;
					 
		} else {
			discarded1++;
			// write to output file
			if (debug) write_to_fastq(discarded_out1, info1, read1, extra1, new_quality1);

			// write to junk output file
			if (debug) write_to_fastq(junk_discarded_out1, info1, trimmed_junk1, extra1, junk_quality1);

			read1_good = false;
		}

		// read2
		if (read2.length() >= minimum_read_length) {
			// write to output file
			if (debug) write_to_fastq(out2, info2, read2, extra2, new_quality2);

			// write to junk output file
			if (debug) write_to_fastq(junk_out2, info2, trimmed_junk2, extra2, junk_quality2);

			read2_good = true;
					 
		} else {
			discarded2++;
			// write to output file
			if (debug) write_to_fastq(discarded_out2, info2, read2, extra2, new_quality2);

			// write to junk output file
			if (debug) write_to_fastq(junk_discarded_out2, info2, trimmed_junk2, extra2, junk_quality2);

			read2_good = false;
		}

		bool concat_good = false;

		// try to concatenate
		if (read1_good && read2_good) {

			string concatenated, concat_quality;

			// concatenate reads
			concatenated = concat_reads(read1, read2, quality1, new_quality2, concat_quality);

			if (concatenated.length() > 0) {
				num_concat++;

				// check the quality of the concatenated read
				if (quality_check(concat_quality)) {
					num_final++;

					string new_info = info1 + string("_concatenated");
					write_to_fastq(concat_out, new_info, concatenated, extra2, concat_quality);

					concat_good = true;
				}
			}
		}

		if (!concat_good) {
			// output the two reads separately
			if (read1_good && quality_check(quality1)) {

				num_non_concat1++;
				write_to_fastq(non_concat_out1, info1, read1, extra2, quality1);

			} else if (read2_good && quality_check(new_quality2)) {

				num_non_concat2++;
				write_to_fastq(non_concat_out2, info2, read2, extra2, new_quality2);

			}
		}
	
	}


	// end timer
	clock_t end_time = clock();

	// find elapsed time
	double elapsed_time = (end_time - start_time) / (double) CLOCKS_PER_SEC;


	in1.close();
	in2.close();

	out1.close();
	junk_out1.close();
	discarded_out1.close();
	junk_discarded_out1.close();
	out2.close();
	junk_out2.close();
	discarded_out2.close();
	junk_discarded_out2.close();


	if (!debug) {
		remove(out_filename1.c_str());
		remove(junk_filename1.c_str());
		remove(discarded_filename1.c_str());
		remove(junk_discarded_filename1.c_str());
		remove(out_filename2.c_str());
		remove(junk_filename2.c_str());
		remove(discarded_filename2.c_str());
		remove(junk_discarded_filename2.c_str());
	}

	concat_out.close();
	non_concat_out1.close();
	non_concat_out2.close();

	
	if (rank == 0) {

		/* merge files */
		if (debug) {
			string out_filename1 = get_file_name(string("trimmed"), 1, string(outfile));
			string junk_filename1 = get_file_name(string("junk"), 1, string(outfile));
			string discarded_filename1 = get_file_name(string("discarded"), 1, string(outfile));
			string junk_discarded_filename1 = get_file_name(string("junk_discarded"), 1, string(outfile));

			string out_filename2 = get_file_name(string("trimmed"), 2, string(outfile));
			string junk_filename2 = get_file_name(string("junk"), 2, string(outfile));
			string discarded_filename2 = get_file_name(string("discarded"), 2, string(outfile));
			string junk_discarded_filename2 = get_file_name(string("junk_discarded"), 2, string(outfile));

			merge_worker_files(out_filename1, size);
			merge_worker_files(junk_filename1, size);
			merge_worker_files(discarded_filename1, size);
			merge_worker_files(junk_discarded_filename1, size);
			merge_worker_files(out_filename2, size);
			merge_worker_files(junk_filename2, size);
			merge_worker_files(discarded_filename2, size);
			merge_worker_files(junk_discarded_filename2, size);
		}

		string non_concat_filename1 = get_file_name(string("non_concat"), 1, string(outfile));
		string non_concat_filename2 = get_file_name(string("non_concat"), 2, string(outfile));
		string concat_filename = string("concat_") + string(outfile);

		merge_worker_files(non_concat_filename1, size);
		merge_worker_files(non_concat_filename2, size);
		merge_worker_files(concat_filename, size);

		
		// get the total numbers for diagnostics
		for (int i = 1; i < size; i++) {
			discarded1 = gather_number(discarded1, i, NUM_DISCARDED1_TAG);
			discarded2 = gather_number(discarded2, i, NUM_DISCARDED2_TAG);
			num_concat = gather_number(num_concat, i, NUM_CONCAT_TAG);
			num_final = gather_number(num_final, i, NUM_FINAL_TAG);
			num_non_concat1 = gather_number(num_non_concat1, i, NUM_NON_CONCAT1_TAG);
			num_non_concat2 = gather_number(num_non_concat2, i, NUM_NON_CONCAT2_TAG);
		}


		// end timer
		clock_t end_time = clock();

		// find elapsed time
		double elapsed_time = (end_time - start_time) / (double) CLOCKS_PER_SEC;

		/* print */
		print_diagnostics(elapsed_time, num_reads1, num_reads2, discarded1, discarded2, num_concat, num_final, num_non_concat1, num_non_concat2);

	} else {
		MPI_Request request;
		MPI_Isend(&discarded1, 1, MPI_INT, 0, NUM_DISCARDED1_TAG, MPI_COMM_WORLD, &request);
		MPI_Isend(&discarded2, 1, MPI_INT, 0, NUM_DISCARDED2_TAG, MPI_COMM_WORLD, &request);
		MPI_Isend(&num_concat, 1, MPI_INT, 0, NUM_CONCAT_TAG, MPI_COMM_WORLD, &request);
		MPI_Isend(&num_final, 1, MPI_INT, 0, NUM_FINAL_TAG, MPI_COMM_WORLD, &request);
		MPI_Isend(&num_non_concat1, 1, MPI_INT, 0, NUM_NON_CONCAT1_TAG, MPI_COMM_WORLD, &request);
		MPI_Isend(&num_non_concat2, 1, MPI_INT, 0, NUM_NON_CONCAT2_TAG, MPI_COMM_WORLD, &request);
	}
}
