/************************************************************
 ****
 ****	This portion reads in a fastq file and fa file
 ****	and trims the adapter sequences from the reads
 ****	in the fastq file. It outputs an output file
 ****	with the trimmed reads along with their quality
 ****	control strings.
 ****
 ****	This version is for runtime testing
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
#include "MPI_readFastq2.hpp"

using namespace std;

// tags for MPI messages
// const int DATA_TAG = 0;
const int DISCARD_TAG = 0;
const int INFO_TAG = 1;
const int SEQ_TAG = 2;
const int EXTRA_TAG = 3;
const int QUAL_TAG = 4;
const int JUNK_SEQ_TAG = 5;
const int JUNK_QUAL_TAG = 6;
const int END_TAG = 10;
const int KEY_TAG = 11;
const int HASH_INFO_TAG = 12;
const int HASH_SEQ_TAG = 13;
const int HASH_QUAL_TAG = 14;
const int NUM_CONCAT_TAG = 15;
const int NUM_FINAL_TAG = 16;


// Hash table
hash_table *myMap;

int NUM_LINES;
int num_reads1, num_reads2, discarded1, discarded2, num_concat, num_final;

// receive char* from MPI and convert it into string
string MPI_receive_string (int source, int tag) {

	int nlength;	// variable for getting string length from MPI
	char* buffer;	// buffer for the received message

	MPI_Status status;
	MPI_Probe(source, tag, MPI_COMM_WORLD, &status);
	MPI_Get_count(&status, MPI_CHAR, &nlength);		// get the size of message

	// allocate buffer
	buffer = (char*) malloc(nlength);

	// receive message
	MPI_Recv(buffer, nlength, MPI_CHAR, source, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	string ret(buffer);

	// free buffer
	free(buffer);

	return ret;
}


// single-threaded
// trim one read file and put the data into the hash table; single-threaded
void trim_file (char* infile, char* outfile, int file_num, vector<string> adapters, bool debug) {

	discarded1 = 0;
	num_reads1 = NUM_LINES;

	ifstream in(infile);
	// open and read from gzipped file
	// ifstream file(infile, ios_base::in | ios_base::binary);
	// boost::iostreams::filtering_istream in;
	// in.push(boost::iostreams::gzip_decompressor());
	// in.push(file);


	// diagnosic files
	string out_filename = get_file_name(string("trimmed"), file_num, string(outfile));
	ofstream out(out_filename);
	string junk_filename = get_file_name(string("junk"), file_num, string(outfile));
	ofstream junk_out(junk_filename);
	string discarded_filename = get_file_name(string("discarded"), file_num, string(outfile));
	ofstream discarded_out(discarded_filename);
	string junk_discarded_filename = get_file_name(string("junk_discarded"), file_num, string(outfile));
	ofstream junk_discarded_out(junk_discarded_filename);


	string info;
	for(int i = 0; i < NUM_LINES; i++) {
	// while (getline(in, info)) {
		getline(in, info);

		if (i+1 % LINEBLOCKS == 0) {
			// cout << "\tParsing read " << i+1 << "..." << endl;
		}

		string sequence, extra, quality;
		string read, trimmed_junk, junk_quality, new_quality;

		// read lines from the input file
		// getline(in, info);		// get fastq read info
		getline(in, sequence);	// get read
		getline(in, extra);		// get the third line
		getline(in, quality);	// get the quality scores 


		// change quality for edge Ts
		strip_t(sequence, quality);

		// trim
		read = trim_read (sequence, quality, new_quality, trimmed_junk, junk_quality, adapters);

		/* write debugging files to outputs */
		if (read.length() >= minimum_read_length) {

			// add into hash table
			string key = myMap->get_key(info);
			myMap->add(key, info, read, new_quality);

			// write to output file
			if (debug) write_to_fastq(out, info, read, extra, new_quality);

			// write to junk output file
			if (debug) write_to_fastq(junk_out, info, trimmed_junk, extra, junk_quality);

		} else {
			discarded1++;

			// write to output file
			if (debug) write_to_fastq(discarded_out, info, read, extra, new_quality);

			// write to junk output file
			if (debug) write_to_fastq(junk_discarded_out, info, trimmed_junk, extra, junk_quality);

		}
	}

	// file.close();
	in.close();

	out.close();
	junk_out.close();
	discarded_out.close();
	junk_discarded_out.close();

	if (!debug) {
		remove(out_filename.c_str());
		remove(junk_filename.c_str());
		remove(discarded_filename.c_str());
		remove(junk_discarded_filename.c_str());
	}
}


// using master/slave parallel processing
// trim one read file and put the data into the hash table; multi-threaded
void MPI_trim_file (char *infile, char *outfile, int file_num, int rank, int size, vector<string> adapters, bool debug) {

	num_reads1 = 0;
	discarded1 = 0;

	int counter = 0;

	if (rank == 0) {

		// diagnosic files
		string out_filename = get_file_name(string("trimmed"), file_num, string(outfile));
		ofstream out(out_filename);
		string junk_filename = get_file_name(string("junk"), file_num, string(outfile));
		ofstream junk_out(junk_filename);
		string discarded_filename = get_file_name(string("discarded"), file_num, string(outfile));
		ofstream discarded_out(discarded_filename);
		string junk_discarded_filename = get_file_name(string("junk_discarded"), file_num, string(outfile));
		ofstream junk_discarded_out(junk_discarded_filename);

		bool continuing = true;
		vector<bool> rank_cont(size, true);

		MPI_Status status;

		while (continuing) {

			// check if worker is done
			int cont;
			MPI_Recv(&cont, 1, MPI_INT, MPI_ANY_SOURCE, END_TAG, MPI_COMM_WORLD, &status);

			// get source of the message
			int curr_source = status.MPI_SOURCE;

			if (!cont) {

				// set worker as finished if end signal is received
				rank_cont[curr_source] = false;

			} else {

				if (num_reads1++ % LINEBLOCKS == 0) {
					// cout << "\tParsing read " << num_reads1 << "..." << endl;
				}

				int discarding;
				char* c_info, c_seq, c_extra, c_qual, c_junk_seq, c_junk_qual;
				string info, read, extra, new_quality, trimmed_junk, junk_quality;

				/* receive messages */
				// find whether if the read is being discarded or not
				MPI_Recv(&discarding, 1, MPI_INT, curr_source, DISCARD_TAG, MPI_COMM_WORLD, &status);

				// get string messages
				info = MPI_receive_string(curr_source, INFO_TAG);
				read = MPI_receive_string(curr_source, SEQ_TAG);
				extra = MPI_receive_string(curr_source, EXTRA_TAG);
				new_quality = MPI_receive_string(curr_source, QUAL_TAG);
				trimmed_junk = MPI_receive_string(curr_source, JUNK_SEQ_TAG);
				junk_quality = MPI_receive_string(curr_source, JUNK_QUAL_TAG);


				/* write debugging files to outputs */
				if (!discarding) {

					// add into hash table
					string key = myMap->get_key(info);
					myMap->add(key, info, read, new_quality);

					// write to output file
					if (debug) write_to_fastq(out, info, read, extra, new_quality);

					// write to junk output file
					if (debug) write_to_fastq(junk_out, info, trimmed_junk, extra, junk_quality);

				} else {

					discarded1++;

					// write to output file
					if (debug) write_to_fastq(discarded_out, info, read, extra, new_quality);

					// write to junk output file
					if (debug) write_to_fastq(junk_discarded_out, info, trimmed_junk, extra, junk_quality);

				}
			}

			// check if all workers are done
			int test;
			for (test = 1; test < size; test++) {
				if (rank_cont[test]) break;		// don't need to check anymore if at least one worker is continuing
			}
			if (test == size) continuing = false;
		}

		out.close();
		junk_out.close();
		discarded_out.close();
		junk_discarded_out.close();

		if (!debug) {
			remove(out_filename.c_str());
			remove(junk_filename.c_str());
			remove(discarded_filename.c_str());
			remove(junk_discarded_filename.c_str());
		}

	} else {

		ifstream in(infile);
		// open and read from gzipped file
		// ifstream file(infile, ios_base::in | ios_base::binary);
		// boost::iostreams::filtering_istream in;
		// in.push(boost::iostreams::gzip_decompressor());
		// in.push(file);


		const int dest = 0;		// always send messages to the master
		int is_continuing = 1;	// signal for the master that this worker is continuing

		for (int i = 0; i < NUM_LINES; i++) {
		string info;
		// while (getline(in, info)) {
		getline(in, info);

			string sequence, extra, quality;
			string read, trimmed_junk, junk_quality, new_quality;

			// read lines from the input file
			// getline(in, info); // get fastq read info
			getline(in, sequence); // get read
			getline(in, extra); // N/A
			getline(in, quality); // get the quality scores 

			if (counter++ % (size-1) == (rank-1) ) {	// only the workers

				//change quality for edge Ts
				strip_t(sequence, quality);

				// trim
				read = trim_read (sequence, quality, new_quality, trimmed_junk, junk_quality, adapters);
				
				int is_discarding;

				if (read.length() >= minimum_read_length) {
					is_discarding = 0;
				} else {
					is_discarding = 1;
				}

				/* send messages */
				// tell the master to continue
				MPI_Send(&is_continuing, 1, MPI_INT, dest, END_TAG, MPI_COMM_WORLD);

				// send data
				MPI_Send(&is_discarding, 1, MPI_INT, dest, DISCARD_TAG, MPI_COMM_WORLD);
				MPI_Send(info.c_str(), info.length()+1, MPI_CHAR, dest, INFO_TAG, MPI_COMM_WORLD);
				MPI_Send(read.c_str(), read.length()+1, MPI_CHAR, dest, SEQ_TAG, MPI_COMM_WORLD);
				MPI_Send(extra.c_str(), extra.length()+1, MPI_CHAR, dest, EXTRA_TAG, MPI_COMM_WORLD);
				MPI_Send(new_quality.c_str(), new_quality.length()+1, MPI_CHAR, dest, QUAL_TAG, MPI_COMM_WORLD);
				MPI_Send(trimmed_junk.c_str(), trimmed_junk.length()+1, MPI_CHAR, dest, JUNK_SEQ_TAG, MPI_COMM_WORLD);
				MPI_Send(junk_quality.c_str(), junk_quality.length()+1, MPI_CHAR, dest, JUNK_QUAL_TAG, MPI_COMM_WORLD);
			}
		}

		// file.close();
		in.close();

		// MPI_Barrier(MPI_COMM_WORLD);

		// signal to master that worker has finished
		is_continuing = 0;
		MPI_Send(&is_continuing, 1, MPI_INT, dest, END_TAG, MPI_COMM_WORLD);
	}

}


// trim reads in file2 and match them with reads from file1 at the same time
// using master/slave parallel processing
void MPI_trim_and_match (char *infile, char *outfile, int file_num, int rank, int size, vector<string> adapters, bool debug) {

	num_reads2 = 0;
	discarded2 = 0;
	num_concat = 0;
	num_final = 0;

	int counter = 0;

	// file for final output file
	string temp = string("final_") + string(outfile);
	stringstream ss;
	ss << rank << "_" << temp;
	string final_out_filename = ss.str();


	if (rank == 0) {

		// diagnosic files
		// each thread would have its own output file
		string out_filename = get_file_name(string("trimmed"), file_num, string(outfile));
		ofstream out(out_filename.c_str());
		string junk_filename = get_file_name(string("junk"), file_num, string(outfile));
		ofstream junk_out(junk_filename.c_str());
		string discarded_filename = get_file_name(string("discarded"), file_num, string(outfile));
		ofstream discarded_out(discarded_filename.c_str());
		string junk_discarded_filename = get_file_name(string("junk_discarded"), file_num, string(outfile));
		ofstream junk_discarded_out(junk_discarded_filename.c_str());


		MPI_Status status;


		// flag for the master to continue
		bool continuing = true;
		vector<bool> rank_cont(size, true);
		while (continuing) {

			// check if worker is done
			int cont;
			MPI_Recv(&cont, 1, MPI_INT, MPI_ANY_SOURCE, END_TAG, MPI_COMM_WORLD, &status);

			// get source of the message
			int curr_source = status.MPI_SOURCE;

			if (!cont) {

				// set worker as finished if end signal is received
				rank_cont[curr_source] = false;

				// get numbers for diagnostics
				int temp_concat, temp_final;
				MPI_Recv(&temp_concat, 1, MPI_INT, curr_source, NUM_CONCAT_TAG, MPI_COMM_WORLD, &status);
				MPI_Recv(&temp_final, 1, MPI_INT, curr_source, NUM_FINAL_TAG, MPI_COMM_WORLD, &status);
				num_concat += temp_concat;
				num_final += temp_final;

			} else {

				if (num_reads2++ % LINEBLOCKS == 0) {
					// cout << "\tParsing read " << num_reads2 << "..." << endl;
				}

				int discarding;
				string info, read, extra, new_quality, trimmed_junk, junk_quality;

				/* receive messages */
				// find whether if the read is being discarded or not
				MPI_Recv(&discarding, 1, MPI_INT, curr_source, DISCARD_TAG, MPI_COMM_WORLD, &status);

				// get string messages
				info = MPI_receive_string(curr_source, INFO_TAG);
				read = MPI_receive_string(curr_source, SEQ_TAG);
				extra = MPI_receive_string(curr_source, EXTRA_TAG);
				new_quality = MPI_receive_string(curr_source, QUAL_TAG);
				trimmed_junk = MPI_receive_string(curr_source, JUNK_SEQ_TAG);
				junk_quality = MPI_receive_string(curr_source, JUNK_QUAL_TAG);


				if (discarding == 0) {

					int contains_key = 1;

					/* move onto the next step (concatenation) */
					// find key
					string key = get_key(info);
					bool found = myMap->has_key(key);


					if (!found) {

						// tell the worker to skip concatenation
						contains_key = 0;
						MPI_Send(&contains_key, 1, MPI_INT, curr_source, KEY_TAG, MPI_COMM_WORLD);

					} else {

						// get value from hash table
						string t_info = myMap->get_info(key);
						string t_seq = myMap->get_seq(key);
						string t_qual = myMap->get_qual(key);

						// erase read1 key from hash table after it's found and used
						myMap->erase(key);

						// tell the worker to continue to the next step
						MPI_Send(&contains_key, 1, MPI_INT, curr_source, KEY_TAG, MPI_COMM_WORLD);
						MPI_Send(t_info.c_str(), t_info.length()+1, MPI_CHAR, curr_source, HASH_INFO_TAG, MPI_COMM_WORLD);
						MPI_Send(t_seq.c_str(), t_seq.length()+1, MPI_CHAR, curr_source, HASH_SEQ_TAG, MPI_COMM_WORLD);
						MPI_Send(t_qual.c_str(), t_qual.length()+1, MPI_CHAR, curr_source, HASH_QUAL_TAG, MPI_COMM_WORLD);
					}

					// write to output file
					if (debug) write_to_fastq(out, info, read, extra, new_quality);

					// write to junk output file
					if (debug) write_to_fastq(junk_out, info, trimmed_junk, extra, junk_quality);

				} else {

					discarded2++;
					// write to output file
					if (debug) write_to_fastq(discarded_out, info, read, extra, new_quality);

					// write to junk output file
					if (debug) write_to_fastq(junk_discarded_out, info, trimmed_junk, extra, junk_quality);
				}

			}

			// check if all workers are done
			int test;
			for (test = 1; test < size; test++) {
				if (rank_cont[test]) break;		// don't need to check anymore if at least one worker is continuing
			}
			if (test == size) continuing = false;
		
		}


		out.close();
		junk_out.close();
		discarded_out.close();
		junk_discarded_out.close();

		
		if (!debug) {
			remove(out_filename.c_str());
			remove(junk_filename.c_str());
			remove(discarded_filename.c_str());
			remove(junk_discarded_filename.c_str());
		}
		

		// check if hash table is empty
		bool empty_hash_table = myMap->is_empty();
		if (!empty_hash_table) {

			ofstream final_out(final_out_filename.c_str());

			// if it's not, put the rest of the reads from file1 into the output file
			for (auto it : myMap->myMap) {
				string curr_pos = it.first;
				string curr_info = it.second.info;
				string curr_seq = it.second.sequence;
				string curr_qual = it.second.quality;

				// quality check
				if (quality_check(curr_qual)) {
					// write to file
					write_to_fastq(final_out, curr_info, curr_seq, string("+"), curr_qual);
				}
			}

			// clear hash table
			myMap->clear();

			final_out.close();
		}

		// cout << "Merging output files..." << endl;
		merge_files(temp, size);

	} else {

		const int dest = 0;	// the destination of the workers' messages is alway going to be the master (rank 0)
		int is_continuing = 1;

		ofstream final_out(final_out_filename.c_str());

		// go through input file2
		ifstream in(infile);
		// open and read from gzipped file
		// ifstream file(infile, ios_base::in | ios_base::binary);
		// boost::iostreams::filtering_istream in;
		// in.push(boost::iostreams::gzip_decompressor());
		// in.push(file);


		for (int i = 0; i < NUM_LINES; i++) {
		string info;
		// while (getline(in, info)) {
		getline(in, info);
			
			string sequence, extra, quality;
			string read, new_quality, trimmed_junk, junk_quality;
			
			// read lines from in1
			getline(in, sequence);
			getline(in, extra);
			getline(in, quality);

			// process only if line is assigned
			if (counter++ % (size-1) == (rank-1) ) {

				//get key for hash table
				string key = get_key(info);

				// change quality for edge Ts
				strip_t(sequence, quality);

				// trim
				read = trim_read (sequence, quality, new_quality, trimmed_junk, junk_quality, adapters);

				int is_discarding;

				// if the trimmed read is valid
				if (read.length() >= minimum_read_length) {
					is_discarding = 0;
				} else {
					is_discarding = 1;
				}

				// send message to the master
				MPI_Send(&is_continuing, 1, MPI_INT, dest, END_TAG, MPI_COMM_WORLD);

				// send data
				MPI_Send(&is_discarding, 1, MPI_INT, dest, DISCARD_TAG, MPI_COMM_WORLD);
				MPI_Send(info.c_str(), info.length()+1, MPI_CHAR, dest, INFO_TAG, MPI_COMM_WORLD);
				MPI_Send(read.c_str(), read.length()+1, MPI_CHAR, dest, SEQ_TAG, MPI_COMM_WORLD);
				MPI_Send(extra.c_str(), extra.length()+1, MPI_CHAR, dest, EXTRA_TAG, MPI_COMM_WORLD);
				MPI_Send(new_quality.c_str(), new_quality.length()+1, MPI_CHAR, dest, QUAL_TAG, MPI_COMM_WORLD);
				MPI_Send(trimmed_junk.c_str(), trimmed_junk.length()+1, MPI_CHAR, dest, JUNK_SEQ_TAG, MPI_COMM_WORLD);
				MPI_Send(junk_quality.c_str(), junk_quality.length()+1, MPI_CHAR, dest, JUNK_QUAL_TAG, MPI_COMM_WORLD);
				

				// continue if read is not discarded
				if (is_discarding == 0) {

					// receive order from master
					int contains_key;
					MPI_Recv(&contains_key, 1, MPI_INT, dest, KEY_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

					// try to concatenate the two reads together if key is found
					if (contains_key == 1) {

						string concat_read, concat_quality;
						string info1, read1, quality1;

						// get messages from master
						info1 = MPI_receive_string(dest, HASH_INFO_TAG);
						read1 = MPI_receive_string(dest, HASH_SEQ_TAG);
						quality1 = MPI_receive_string(dest, HASH_QUAL_TAG);

						// try to concatenate
						concat_read = concat_reads(read1, read, quality1, new_quality, concat_quality);


						// if concatenation is possible
						if (concat_read.length() > 0) {

							num_concat++;

							// output concatenated read
							if (quality_check(concat_quality)) {
								num_final++;

								string new_info = info + string("concatenated_");
								write_to_fastq(final_out, new_info, concat_read, extra, concat_quality);
							}

						} else {	// otherwise

							// output the two reads separately
							// check quality1
							if (quality_check(quality1)) {
								write_to_fastq(final_out, info1, read1, extra, quality1);
							}

							// check quality2
							if (quality_check(new_quality)) {
								write_to_fastq(final_out, info, read, extra, new_quality);
							}

						}
					} else {

						// output read2 on its own
						if (quality_check(new_quality)) {
							write_to_fastq(final_out, info, read, extra, new_quality);
						}

					}
				}
			}
			
			// counter++;
		}

		// file.close();
		in.close();

		final_out.close();

		// signal to master that the worker has finished
		is_continuing = 0;
		MPI_Send(&is_continuing, 1, MPI_INT, dest, END_TAG, MPI_COMM_WORLD);
		MPI_Send(&num_concat, 1, MPI_INT, dest, NUM_CONCAT_TAG, MPI_COMM_WORLD);
		MPI_Send(&num_final, 1, MPI_INT, dest, NUM_FINAL_TAG, MPI_COMM_WORLD);
	}
}


// single threaded read processing function
double process_reads (char* infile1, char* infile2, char* outfile, int file_num, vector<string> adapters, int lines, bool debug) {

	NUM_LINES = lines;
	num_reads2 = NUM_LINES;
	discarded1 = 0; discarded2 = 0; num_concat = 0; num_final = 0;

	// make a hash table
	myMap = new hash_table();

	//start clock
	clock_t start_time = clock();

	// cout << "Reading file1..." << endl;

	// get hash table for read file1
	trim_file(infile1, outfile, 1, adapters, debug);

	// files for infile2
	string out_filename = get_file_name(string("trimmed"), file_num, string(outfile));
	ofstream out2(out_filename);
	string junk_filename = get_file_name(string("junk"), file_num, string(outfile));
	ofstream junk_out2(junk_filename);
	string discarded_filename = get_file_name(string("discarded"), file_num, string(outfile));
	ofstream discarded_out2(discarded_filename);
	string junk_discarded_filename = get_file_name(string("junk_discarded"), file_num, string(outfile));
	ofstream junk_discarded_out2(junk_discarded_filename);

	// cout << "Reading file2..." << endl;

	// file for concatenated reads
	string final_out_filename = string("final_") + string(outfile);
	ofstream final_out(final_out_filename.c_str());

	ifstream in2(infile2);
	// open and read from gzipped file
	// ifstream file2(infile2, ios_base::in | ios_base::binary);
	// boost::iostreams::filtering_istream in2;
	// in2.push(boost::iostreams::gzip_decompressor());
	// in2.push(file2);


	string info2;
	for (int i = 0; i < NUM_LINES; i++) {
	// while (getline(in2, info2)) {
		getline(in2, info2);

		if (i+1 % LINEBLOCKS == 0) {
			// cout << "\tParsing read " << i+1 << "..." << endl;
		}
		
		string sequence2, extra2, quality2;
		string read2, trimmed_junk2, junk_quality2, new_quality2;

		// file 2
		getline(in2, sequence2); // get read 2
		getline(in2, extra2); // N/A
		getline(in2, quality2); // get the quality score 2

		//strip T's for file 1
		strip_t(sequence2, quality2);
		
		// trim read 2
		read2 = trim_read (sequence2, quality2, new_quality2, trimmed_junk2, junk_quality2, adapters);

		/* write debugging files to outputs */
		if (read2.length() >= minimum_read_length) {
			// write to output file
			if (debug) write_to_fastq(out2, info2, read2, extra2, new_quality2);

			// write to junk output file
			if (debug) write_to_fastq(junk_out2, info2, trimmed_junk2, extra2, junk_quality2);
					 
		} else {
			discarded2++;
			// write to output file
			if (debug) write_to_fastq(discarded_out2, info2, read2, extra2, new_quality2);

			// write to junk output file
			if (debug) write_to_fastq(junk_discarded_out2, info2, trimmed_junk2, extra2, junk_quality2);
		}

		if (read2.length() >= minimum_read_length) {

			// get the key for read2
			string key = myMap->get_key(info2);
			bool found = myMap->has_key(key);

			// check if the hash table has said key
			if (!found) {

				if (quality_check(new_quality2)) {

					// output read2 on its own
					write_to_fastq(final_out, info2, read2, extra2, new_quality2);
				}

			} else {

				string info1, read1, quality1, concatenated, concat_quality;

				info1 = myMap->get_info(key);
				read1 = myMap->get_seq(key);
				quality1 = myMap->get_qual(key);

				myMap->erase(key);

				// concatenate reads
				concatenated = concat_reads(read1, read2, quality1, new_quality2, concat_quality);

				if (concatenated.length() > 0) {
					num_concat++;

					// concat_out << info1 << endl
					// 		   << concatenated << endl
					// 		   << extra1 << endl
					// 		   << new_quality << endl;

					// check the quality of the concatenated read
					if (quality_check(concat_quality)) {
						num_final++;

						string new_info = info1 + string("concatenated_");
						write_to_fastq(final_out, new_info, concatenated, extra2, concat_quality);
					}

				} else {

					// output the two reads separately
					if (quality_check(quality1)) {
						write_to_fastq(final_out, info1, read1, extra2, quality1);
					}

					if (quality_check(new_quality2)) {
						write_to_fastq(final_out, info2, read2, extra2, new_quality2);
					}
				}
			}
		}
	}


	// check if hash table is empty
	bool empty_hash_table = myMap->is_empty();
	if (!empty_hash_table) {

		// if it's not, put the rest of the reads from file1 into the output file
		for (auto it : myMap->myMap) {
			string curr_pos = it.first;
			string curr_info = it.second.info;
			string curr_seq = it.second.sequence;
			string curr_qual = it.second.quality;

			// quality check
			if (quality_check(curr_qual)) {
				write_to_fastq(final_out, curr_info, curr_seq, string("+"), curr_qual);
			}
		}

		// clear hash table
		myMap->clear();
	}


	// end timer
	clock_t end_time = clock();

	// find elapsed time
	double elapsed_time = (end_time - start_time) / (double) CLOCKS_PER_SEC;

	// print_diagnostics(elapsed_time, num_reads1, j, discarded1, discarded2, num_concat, num_final);

	// file2.close();
	in2.close();

	out2.close();
	junk_out2.close();
	discarded_out2.close();
	junk_discarded_out2.close();


	if (!debug) {
		remove(out_filename.c_str());
		remove(junk_filename.c_str());
		remove(discarded_filename.c_str());
		remove(junk_discarded_filename.c_str());
	}


	final_out.close();

	// delete hash table
	delete myMap;

	return elapsed_time;
}


double MPI_process_reads (char* infile1, char* infile2, char* outfile, int rank, int size, vector<string> adapters, int lines, bool debug) {

	NUM_LINES = lines;

	discarded1 = 0; discarded2 = 0; num_reads1 = 0; num_reads2 = 0; num_concat = 0; num_final = 0;

	// make a hash table
	myMap = new hash_table();

	//start clock
	// clock_t start_time = clock();
	double start_time = MPI_Wtime();

	if (rank == 0) {
		// cout << "Reading file1..." << endl;
	}

	// trim read 1
	MPI_trim_file(infile1, outfile, 1, rank, size, adapters, debug);	// hash table would only exist in the master process (rank 0)

	// wait for every process to finish
	MPI_Barrier(MPI_COMM_WORLD);

	if (rank == 0) {
		// cout << "Reading file2..." << endl;
	}

	// trim read 2 and merge read 1 and 2 together
	MPI_trim_and_match(infile2, outfile, 2, rank, size, adapters, debug);

	// wait for every process to finish
	MPI_Barrier(MPI_COMM_WORLD);

	// delete hash table
	delete myMap;

	// end timer
	// clock_t end_time = clock();
	double end_time = MPI_Wtime();


	if (rank == 0) {
		// print_diagnostics(end_time - start_time, num_reads1, num_reads2, discarded1, discarded2, num_concat, num_final);
	}

	return end_time - start_time;
}

