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
#include <unordered_map>
#include <mpi.h>
#include "SmithWaterman.hpp"
#include "sys/types.h"
#include "sys/sysinfo.h"
#include "sys/times.h"
#include "sys/vtimes.h"

using namespace std;

#define LINEBLOCKS 500

// defaults
#define DEFAULT_MIN_CUT_LENGTH 11
#define DEFAULT_MIN_MATCH_LENGTH 7

// cutoff values
int minimum_read_length = DEFAULT_MIN_CUT_LENGTH;
int minimum_match_length = DEFAULT_MIN_MATCH_LENGTH;

// universal sequence for tru seq adapters
const string TRU_SEQ_ADAPTER = "GATCGGAAGAGCACACGTCTGAACTCCAGTCAC";

class SmithWaterman;

// containers
vector<vector<string> > file_line;
vector< vector<string> > trimmed_seq;
vector<string> adapters;

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
const int CONCAT_TAG = 15;
const int CONCAT_INFO_TAG = 16;
const int CONCAT_SEQ_TAG = 17;
const int CONCAT_QUAL_TAG = 18;
const int CONCAT_END_TAG = 20;




struct sysinfo memInfo;
static unsigned long long lastTotalUser, lastTotalUserLow, lastTotalSys, lastTotalIdle;
static clock_t lastCPU, lastSysCPU, lastUserCPU;
static int numProcessors;

// struct for holding info for a read
struct DNAread {
	string sequence;
	string quality;
	string info;
};

// Hash table
unordered_map<string, DNAread> myMap;


// change quality score to 2 ('#') if the leading and/or last 2 bases are 'T's
bool strip_t (string &sequence, string &quality) {
	bool has_t = false;
	int seq_length = sequence.length();
	if (sequence[0] == 'T') {
		quality[0] = '#';
		has_t = true;
	}
	if (sequence[seq_length-2] == 'T' && sequence[seq_length-1] == 'T') {
		quality[seq_length-2] = quality[seq_length-1] = '#';
		has_t = true;
	}
	return has_t;
}


// finds bases with score less than 15 (more than 5: fails)
bool quality_check (string quality) {
	int num_failed = 0;
	for (int i = 0; i < quality.length(); i++) {
		char cur = quality[i];
		if (int(cur) - 33 < 15) {	// adjust ascii value for base 33 quality strings
			num_failed++;
		}
		if (num_failed > 5) return false;
	}
	return true;
}


string trim_read (string sequence, string quality, string &new_quality, string &junk, string &junk_quality, bool has_adapter_file) {
	string read = sequence;
	int highest_score = 0;

	if (has_adapter_file) {

		for (int i = 0; i < adapters.size(); i++) {

			SmithWaterman* sw = new SmithWaterman(sequence, adapters[i], quality, "", minimum_match_length);

			int curr_high_score = sw->get_highest_score();

			if (highest_score < curr_high_score) {
				highest_score = curr_high_score;
				read = sw->trim_from_ending();
				new_quality = sw->get_quality1();
				junk = sw->get_trimmed();
				junk_quality = sw->get_trimmed_quality();
			}

			delete sw;
		}
	} else {

		SmithWaterman* sw = new SmithWaterman(sequence, TRU_SEQ_ADAPTER, quality, "", minimum_match_length);

		int curr_high_score = sw->get_highest_score();

		if (highest_score < curr_high_score) {
			highest_score = curr_high_score;
			read = sw->trim_from_ending();
			new_quality = sw->get_quality1();
			junk = sw->get_trimmed();
			junk_quality = sw->get_trimmed_quality();
		}

		delete sw;
	}

	return read;
}


string concat_reads (string sequence1, string sequence2, string quality1, string quality2, string &new_quality) {
	string concat = "";
	int highest_score = 0;


	SmithWaterman* sw = new SmithWaterman(sequence1, sequence2, quality1, quality2, minimum_match_length);
	int curr_high_score = sw->get_highest_score();
	if (highest_score < curr_high_score) {
		if (sw->match_reads()) {
			concat = sw->get_matched_string();
			new_quality = sw->get_matched_quality();
		}
	}
	delete sw;

	return concat;
}


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


// get key for hash table
string get_key (string info) {

	// tokenize
	stringstream ss;
	ss.str(info);
	string token;
	for (int i = 0; i < 5; i++) {
		getline(ss, token, ':');
	}
		
	string x;
	getline(ss, x, ':');
	string y;
	getline(ss, y, ' ');
	
	string key = x + ":" + y;

	return key;	
}


// trim one read file and put the data into the hash table; single-threaded
void trim_file (char *infile, char *outfile, int file_num, bool has_adapter_file, int &discarded) {

	discarded = 0;

	ifstream in(infile);

	// diagnosic files
	string out_filename = get_file_name(string("trimmed"), file_num, string(outfile));
	ofstream out(out_filename);
	string junk_filename = get_file_name(string("junk"), file_num, string(outfile));
	ofstream junk_out(junk_filename);
	string discarded_filename = get_file_name(string("discarded"), file_num, string(outfile));
	ofstream discarded_out(discarded_filename);
	string junk_discarded_filename = get_file_name(string("junk_discarded"), file_num, string(outfile));
	ofstream junk_discarded_out(junk_discarded_filename);


	int j = 0;
	string info;
	while (getline(in, info)) {

		if (j++ % LINEBLOCKS == 0) {
			cout << "\tParsing read " << j << "..." << endl;
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
		read = trim_read (sequence, quality, new_quality, trimmed_junk, junk_quality, has_adapter_file);

		/* write debugging files to outputs */
		if (read.length() >= minimum_read_length) {

			// add into hash table
			string key = get_key(info);
			DNAread tempLine;
			tempLine.sequence = read;
			tempLine.quality = new_quality;
			tempLine.info = info;
			myMap.insert({key, tempLine});

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
			// write to output file
			discarded_out << info << endl
						  << read << endl
						  << extra << endl
						  << new_quality << endl;

			// write to junk output file
			junk_discarded_out << info << endl
							   << trimmed_junk << endl
							   << extra << endl
							   << junk_quality << endl;
		}
	}

	out.close();
	junk_out.close();
	discarded_out.close();
	junk_discarded_out.close();
}


// using master/slave parallel processing
// trim one read file and put the data into the hash table; multi-threaded
void MPI_trim_file (char *infile, char *outfile, int file_num, int rank, int size, bool has_adapter_file, int &discarded) {

	discarded = 0;

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

		int j = 0;

		while (continuing) {

			// check if worker is done
			int cont;
			MPI_Recv(&cont, 1, MPI_INT, MPI_ANY_SOURCE, END_TAG, MPI_COMM_WORLD, &status);

			if (j++ % LINEBLOCKS == 0) {
				cout << "\tParsing read " << j << "..." << endl;
			}

			// get source of the message
			int curr_source = status.MPI_SOURCE;

			if (!cont) {

				// set worker as finished if end signal is received
				rank_cont[curr_source] = false;

			} else {

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
					string key = get_key(info);
					DNAread tempLine;
					tempLine.sequence = read;
					tempLine.quality = new_quality;
					tempLine.info = info;
					myMap.insert({key, tempLine});

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

					// write to output file
					discarded_out << info << endl
								  << read << endl
								  << extra << endl
								  << new_quality << endl;

					// write to junk output file
					junk_discarded_out << info << endl
									   << trimmed_junk << endl
									   << extra << endl
									   << junk_quality << endl;

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

	} else {

		ifstream in(infile);

		const int dest = 0;		// always send messages to the master
		int is_continuing = 1;	// signal for the master that this worker is continuing

		string info;
		while (getline(in, info)) {
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
				read = trim_read (sequence, quality, new_quality, trimmed_junk, junk_quality, has_adapter_file);
				
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

		// MPI_Barrier(MPI_COMM_WORLD);

		// signal to master that worker has finished
		is_continuing = 0;
		MPI_Send(&is_continuing, 1, MPI_INT, dest, END_TAG, MPI_COMM_WORLD);
	}

}


/*** WIP -- use master/slave to output file? ***/
void MPI_trim_and_match (char *infile, char *outfile, int file_num, int rank, int size, bool has_adapter_file, int &discarded) {

	discarded = 0;

	int counter = 0;

	// file for concatenated reads
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

		int j = 0;

		// flag for the master to continue
		bool continuing = true;
		vector<bool> rank_cont(size, true);
		while (continuing) {

			// check if worker is done
			int cont;
			MPI_Recv(&cont, 1, MPI_INT, MPI_ANY_SOURCE, END_TAG, MPI_COMM_WORLD, &status);

			if (j++ % LINEBLOCKS == 0) {
				cout << "\tParsing read " << j << "..." << endl;
			}

			// get source of the message
			int curr_source = status.MPI_SOURCE;

			if (!cont) {

				// set worker as finished if end signal is received
				rank_cont[curr_source] = false;

			} else {

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

					int has_key = 1;

					/* move onto the next step (concatenation) */
					// find key
					string key = get_key(info);
					unordered_map<string, DNAread>::const_iterator found = myMap.find(key);


					if (found == myMap.end()) {

						// tell the worker to skip concatenation
						has_key = 0;
						MPI_Send(&has_key, 1, MPI_INT, curr_source, KEY_TAG, MPI_COMM_WORLD);

					} else {

						// get value from hash table
						DNAread value = myMap[key];
						string t_info = value.info;
						string t_seq = value.sequence;
						string t_qual = value.quality;

						// erase read1 key from hash table after it's found and used
						myMap.erase(key);

						// tell the worker to continue to the next step
						MPI_Send(&has_key, 1, MPI_INT, curr_source, KEY_TAG, MPI_COMM_WORLD);
						MPI_Send(t_info.c_str(), t_info.length()+1, MPI_CHAR, curr_source, HASH_INFO_TAG, MPI_COMM_WORLD);
						MPI_Send(t_seq.c_str(), t_seq.length()+1, MPI_CHAR, curr_source, HASH_SEQ_TAG, MPI_COMM_WORLD);
						MPI_Send(t_qual.c_str(), t_qual.length()+1, MPI_CHAR, curr_source, HASH_QUAL_TAG, MPI_COMM_WORLD);
					}

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
					// write to output file
					discarded_out << info << endl
									<< read << endl
									<< extra << endl
									<< new_quality << endl;

					// write to junk output file
					junk_discarded_out << info << endl
										<< trimmed_junk << endl
										<< extra << endl
										<< junk_quality << endl;

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
		

		// check if hash table is empty
		bool empty_hash_table = myMap.empty();
		if (!empty_hash_table) {

			ofstream final_out(final_out_filename.c_str());

			// if it's not, put the rest of the reads from file1 into the output file
			for (auto it : myMap) {
				string curr_pos = it.first;
				string curr_info = it.second.info;
				string curr_seq = it.second.sequence;
				string curr_qual = it.second.quality;

				// quality check
				if (quality_check(curr_qual)) {
					final_out << curr_info << endl
							  << curr_seq << endl
							  << "+" << endl
							  << curr_qual << endl;
				}
			}

			// clear hash table
			myMap.clear();

			final_out.close();
		}

		cout << "Merging output files..." << endl;
		merge_files(temp, size);

	} else {

		const int dest = 0;	// the workers' message destination is alway going to be the master (rank 0)
		int is_continuing = 1;

		ofstream final_out(final_out_filename.c_str());

		// go through input file2
		ifstream in(infile);
		string info;
		while (getline(in, info)) {
			
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
				read = trim_read (sequence, quality, new_quality, trimmed_junk, junk_quality, has_adapter_file);

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
					int has_key;
					MPI_Recv(&has_key, 1, MPI_INT, dest, KEY_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

					// try to concatenate the two reads together if key is found
					if (has_key == 1) {

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

							// output concatenated read
							if (quality_check(concat_quality)) {
								final_out << "concatenated_" << info << endl
										  << concat_read << endl
										  << extra << endl
										  << concat_quality << endl;
							}

						} else {	// otherwise

							// output the two reads separately
							// check quality1
							if (quality_check(quality1)) {
								final_out << info1 << endl
										  << read1 << endl
										  << extra << endl
										  << quality1 << endl;
							}

							// check quality2
							if (quality_check(new_quality)) {
								final_out << info << endl
										  << read << endl
										  << extra << endl
										  << new_quality << endl;
							}

						}
					} else {

						// output read2 on its own
						if (quality_check(new_quality)) {
							final_out << info << endl
									  << read << endl
									  << extra << endl
									  << new_quality << endl;
						}

					}
				}
			}
			
			// counter++;
		}

		final_out.close();

		// signal to master that the worker has finished
		is_continuing = 0;
		MPI_Send(&is_continuing, 1, MPI_INT, dest, END_TAG, MPI_COMM_WORLD);
	}
}


// single threaded read processing function
void process_reads (char* infile1, char* infile2, char* outfile, int file_num, bool has_adapter_file) {

	int discarded1 = 0, discarded2 = 0, num_concat = 0, num_final = 0;

	//start clock
	clock_t start_time = clock();

	cout << "Reading file1..." << endl;

	// get hash table for read file1
	trim_file(infile1, outfile, 1, has_adapter_file, discarded1);


	// files for infile2
	string out_filename = get_file_name(string("trimmed"), file_num, string(outfile));
	ofstream out2(out_filename);
	string junk_filename = get_file_name(string("junk"), file_num, string(outfile));
	ofstream junk_out2(junk_filename);
	string discarded_filename = get_file_name(string("discarded"), file_num, string(outfile));
	ofstream discarded_out2(discarded_filename);
	string junk_discarded_filename = get_file_name(string("junk_discarded"), file_num, string(outfile));
	ofstream junk_discarded_out2(junk_discarded_filename);


	cout << "Reading file2..." << endl;

	// file for concatenated reads
	string final_out_filename = string("final_") + string(outfile);
	ofstream final_out(final_out_filename.c_str());

	ifstream in2(infile2);

	int j = 0;
	string info2;
	while (getline(in2, info2)) {

		if (j++ % LINEBLOCKS == 0) {
			cout << "\tParsing read " << j << "..." << endl;
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
		read2 = trim_read (sequence2, quality2, new_quality2, trimmed_junk2, junk_quality2, has_adapter_file);

		/* write debugging files to outputs */
		if (read2.length() >= minimum_read_length) {
			// write to output file
			out2 << info2 << endl
				<< read2 << endl
				<< extra2 << endl
				<< new_quality2 << endl;

			// write to junk output file
			junk_out2 << info2 << endl
					 << trimmed_junk2 << endl
					 << extra2 << endl
					 << junk_quality2 << endl;
					 
		} else {
			discarded2++;
			// write to output file
			discarded_out2 << info2 << endl
						  << read2 << endl
						  << extra2 << endl
						  << new_quality2 << endl;

			// write to junk output file
			junk_discarded_out2 << info2 << endl
							   << trimmed_junk2 << endl
							   << extra2 << endl
							   << junk_quality2 << endl;
		}

		if (read2.length() >= minimum_read_length) {

			// check if hash table has key
			string key = get_key(info2);
			if (myMap.find(key) == myMap.end()) {

				if (quality_check(read2)) {

					// output read2 on its own
					final_out << info2 << endl
							  << read2 << endl
							  << extra2 << endl
							  << new_quality2 << endl;
				}

			} else {

				bool has_concat = false;
				string info1, read1, quality1, concatenated, concat_quality;

				info1 = myMap[key].info;
				read1 = myMap[key].sequence;
				quality1 = myMap[key].quality;

				myMap.erase(key);

				// concatenate reads
				concatenated = concat_reads(read1, read2, quality1, quality2, concat_quality);

				if (concatenated.length() > 0) {
					num_concat++;
					has_concat = true;

					// concat_out << info1 << endl
					// 		   << concatenated << endl
					// 		   << extra1 << endl
					// 		   << new_quality << endl;

					// check the quality of the concatenated read
					if (quality_check(concat_quality)) {
						num_final++;

						final_out << "concatenated_" << info1 << endl
								  << concatenated << endl
								  << extra2 << endl
								  << concat_quality << endl;
					}

				} else {

					// output the two reads separately
					final_out << info1 << endl
							  << read1 << endl
							  << extra2 << endl
							  << quality1 << endl;

					final_out << info2 << endl
							  << read2 << endl
							  << extra2 << endl
							  << new_quality2 << endl;
				}
			}
		}
	}


	// check if hash table is empty
	bool empty_hash_table = myMap.empty();
	if (!empty_hash_table) {

		// if it's not, put the rest of the reads from file1 into the output file
		for (auto it : myMap) {
			string curr_pos = it.first;
			string curr_info = it.second.info;
			string curr_seq = it.second.sequence;
			string curr_qual = it.second.quality;

			// quality check
			if (quality_check(curr_qual)) {
				final_out << curr_info << endl
						  << curr_seq << endl
						  << "+" << endl
						  << curr_qual << endl;
			}
		}

		// clear hash table
		myMap.clear();
	}


	// end timer
	clock_t end_time = clock();

	// find elapsed time
	double elapsed_time = (end_time - start_time) / (double) CLOCKS_PER_SEC;
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
	cout << "Total number of reads: " << j << endl;
	cout << "Number of discarded reads in file 1 due to being trimmed too short: " << discarded1 << endl;
	cout << "Number of discarded reads in file 2 due to being trimmed too short: " << discarded2 << endl;
	cout << "Number of concatenated reads: " << num_concat << endl;
	cout << "NUmber of concatenated reads that passed quality check: " << num_final << endl;


	in2.close();
	out2.close();
	junk_out2.close();
	discarded_out2.close();
	junk_discarded_out2.close();

	final_out.close();

	return;
}

void MPI_process_reads (char* infile1, char* infile2, char* outfile, int rank, int size, bool has_adapter_file) {

	int discarded1, discarded2;

	//start clock
	// clock_t start_time = clock();
	double start_time = MPI_Wtime();

	if (rank == 0) {
		cout << "Reading file1..." << endl;
	}

	// trim read 1
	MPI_trim_file(infile1, outfile, 1, rank, size, has_adapter_file, discarded1);	// hash table would only exist in the master process (rank 0)
	// trim_file(infile1, outfile, 1, has_adapter_file, discarded1);

	// wait for every process to finish
	MPI_Barrier(MPI_COMM_WORLD);

	if (rank == 0) {
		cout << "Reading file2..." << endl;
	}

	// trim read 2 and merge read 1 and 2 together
	MPI_trim_and_match(infile2, outfile, 2, rank, size, has_adapter_file, discarded2);


	// end timer
	// clock_t end_time = clock();
	double end_time = MPI_Wtime();


	if (rank == 0) {
		// find elapsed time
		double elapsed_time = end_time - start_time;
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


		cout << "Minimum read length after trimming: " << minimum_read_length << endl;
		cout << "Minimum matching length: " << minimum_match_length << endl;
		// cout << "Total number of reads: " << j << endl;
		cout << "Number of discarded reads in file 1 due to being trimmed too short: " << discarded1 << endl;
		cout << "Number of discarded reads in file 2 due to being trimmed too short: " << discarded2 << endl;
		// cout << "Number of concatenated reads: " << num_concat << endl;
		// cout << "NUmber of concatenated reads that passed quality check: " << num_final << endl;
	}
}

void print_argument_directions(char* program) {
	cerr << "Usage: " << program << " <path_to_file1> <path_to_file2> <output_filename> " << endl;
	cerr << "optional: --adapter <path_to_adapter_file> (default: truseq adapters)";
	cerr << "             (*note that using your own adapter file may slow the program down)" << endl;
	cerr << "          --trim <minimum length of trimmed reads> (default: 11)" << endl;
	cerr << "          --match <minimum match length> (default: 7)" << endl;
	cerr << "          --MPI (flag for using MPI)" << endl;
}


int parseLine(char* line){
    // This assumes that a digit will be found and the line ends in " Kb".
    int i = strlen(line);
    const char* p = line;
    while (*p <'0' || *p > '9') p++;
    line[i-3] = '\0';
    i = atoi(p);
    return i;
}

int getValue(){ //Note: this value is in KB!
    FILE* file = fopen("/proc/self/status", "r");
    int result = -1;
    char line[128];

    while (fgets(line, 128, file) != NULL){
        if (strncmp(line, "VmRSS:", 6) == 0){
            result = parseLine(line);
            break;
        }
    }
    fclose(file);
    return result;
}


void initCPUCurrent(){
    FILE* file = fopen("/proc/stat", "r");
    fscanf(file, "cpu %llu %llu %llu %llu", &lastTotalUser, &lastTotalUserLow,
        &lastTotalSys, &lastTotalIdle);
    fclose(file);
}

double CPUCurrent(){
    double percent;
    FILE* file;
    unsigned long long totalUser, totalUserLow, totalSys, totalIdle, total;

    file = fopen("/proc/stat", "r");
    fscanf(file, "cpu %llu %llu %llu %llu", &totalUser, &totalUserLow,
        &totalSys, &totalIdle);
    fclose(file);

    if (totalUser < lastTotalUser || totalUserLow < lastTotalUserLow ||
        totalSys < lastTotalSys || totalIdle < lastTotalIdle){
        //Overflow detection. Just skip this value.
        percent = -1.0;
    }
    else{
        total = (totalUser - lastTotalUser) + (totalUserLow - lastTotalUserLow) +
            (totalSys - lastTotalSys);
        percent = total;
        total += (totalIdle - lastTotalIdle);
        percent /= total;
        percent *= 100;
    }

    lastTotalUser = totalUser;
    lastTotalUserLow = totalUserLow;
    lastTotalSys = totalSys;
    lastTotalIdle = totalIdle;

    cout << percent << endl;

    return percent;
}


void initCPUCurrentProcess(){
    FILE* file;
    struct tms timeSample;
    char line[128];

    lastCPU = times(&timeSample);
    lastSysCPU = timeSample.tms_stime;
    lastUserCPU = timeSample.tms_utime;

    file = fopen("/proc/cpuinfo", "r");
    numProcessors = 0;
    while(fgets(line, 128, file) != NULL){
        if (strncmp(line, "processor", 9) == 0) numProcessors++;
    }
    fclose(file);
}

double	CPUCurrentProcess(){
    struct tms timeSample;
    clock_t now;
    double percent;

    now = times(&timeSample);
    if (now <= lastCPU || timeSample.tms_stime < lastSysCPU ||
        timeSample.tms_utime < lastUserCPU){
        //Overflow detection. Just skip this value.
        percent = -1.0;
    }
    else{
        percent = (timeSample.tms_stime - lastSysCPU) +
            (timeSample.tms_utime - lastUserCPU);
        percent /= (now - lastCPU);
        percent /= numProcessors;
        percent *= 100;
    }
    lastCPU = now;
    lastSysCPU = timeSample.tms_stime;
    lastUserCPU = timeSample.tms_utime;

    cout << percent << endl;

    return percent;
}

void GetMemoryUsage() {
	sysinfo (&memInfo);
	long long totalVirtualMem = memInfo.totalram;
	//Add other values in next statement to avoid int overflow on right hand side...
	totalVirtualMem += memInfo.totalswap;
	totalVirtualMem *= memInfo.mem_unit;
	long long totalPhysMem = memInfo.totalram;
	//Multiply in next statement to avoid int overflow on right hand side...
	totalPhysMem *= memInfo.mem_unit;
	long long physMemUsed = memInfo.totalram - memInfo.freeram;
	//Multiply in next statement to avoid int overflow on right hand side...
	physMemUsed *= memInfo.mem_unit;

	cout << "Total physical memory: " << totalPhysMem << endl;
	cout << "Physical memory currently used" << physMemUsed << endl;
	cout << "CPU currently csed" << initCPUCurrent(); << endl;
	cout << "CPU currently used by current process:" << initCPUCurrentProcess(); << endl;	
}


int main (int argc, char** argv) {

	if (argc < 4) {
		print_argument_directions(argv[0]);
		exit(1); 
	}

	initCPUCurrent();
	initCPUCurrentProcess();

	bool using_adapter_file = false;
	bool use_MPI = false;

	// grab values from arguments
	if (argc > 4) {
		int i = 4;
		while (i < argc) {

			if (string(argv[i]).compare("--adapter") == 0) {
				using_adapter_file = true;
				// get the list of adapters
				ifstream adpt(argv[++i]);
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
			} else if (string(argv[i]).compare("--trim") == 0) {
				minimum_read_length = atoi(argv[++i]);
			} else if (string(argv[i]).compare("--match") == 0) {
				minimum_match_length = atoi(argv[++i]);
			} else if (string(argv[i]).compare("--MPI") == 0) {
				use_MPI = true;
			} else {
				print_argument_directions(argv[0]);
				exit(0);
			}

			i++;
		}
	}

	if (use_MPI) {
		// start MPI
		int size, rank;
		MPI_Init(&argc, &argv);
		MPI_Comm_size(MPI_COMM_WORLD, &size);
		MPI_Comm_rank(MPI_COMM_WORLD, &rank);

		MPI_process_reads(argv[1], argv[2], argv[3], rank, size, using_adapter_file);

		MPI_Finalize();
	} else {
		process_reads(argv[1], argv[2], argv[3], 2, using_adapter_file);
	}

	return 0;
}

