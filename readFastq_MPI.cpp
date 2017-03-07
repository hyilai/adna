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

using namespace std;

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


// struct Coordinate {
// 	int x;
// 	int y;
// };

// Coordinate find_coord (string info) {
// 	stringstream ss;
// 	ss.str(info);
// 	string token;
// 	for (int i = 0; i < 5; i++) {
// 		getline(ss, token, ':');
// 	}

// 	Coordinate pt;
// 	getline(ss, token, ':');
// 	pt.x = atoi(token.c_str());
// 	getline(ss, token, ':');
// 	pt.y = atoi(token.c_str());

// 	return pt;
// }

// bool is_same_coord (Coordinate coord1, Coordinate coord2) {
// 	return (coord1.x == coord2.x && coord1.y == coord2.y);
// }

struct DNAread {
	string sequence;
	string quality;
	string coor;
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

// get key for hash table
string get_key (string info) {

	// tokenize
	stringstream ss;
	ss.str(info);
	string token;
	for (int i = 0; i < 2; i++) {
		getline(ss, token, ':');
	}
		
	string x;
	getline(ss, x, ':');
	string y;
	getline(ss, y, ':');
	
	string key = x + ":" + y;

	return key;	
}


// trim one read file and put the data into the hash table
void trim_file (char *infile, char *outfile, int file_num, bool has_adapter_file) {

	int discarded = 0;

	// diagnosic files
	ifstream in(infile);
	string out_filename = "trimmed_" + file_num + "_" + string(outfile);
	ofstream out(out_filename.c_str());
	string junk_filename = "junk_" + file_num + "_" + string(outfile);
	ofstream junk_out(junk_filename.c_str());
	string discarded_filename = "discarded_" + file_num + "_" + string(outfile);
	ofstream discarded_out(discarded_filename.c_str());
	string junk_discarded_filename = "junk_discarded_" + file_num + "_" + string(outfile);
	ofstream junk_discarded_out(junk_discarded_filename.c_str());

	string info;
	while (getline(in, info)) {
		string sequence, extra, quality;
		string read, trimmed_junk, junk_quality, new_quality;

		// read lines from the input file
		// getline(in, info); // get fastq read info
		getline(in, sequence); // get read
		getline(in, extra); // N/A
		getline(in, quality); // get the quality scores 


		//strip T's for file 1
		strip_t(sequence, quality);

		// trim read 1
		read = trim_read (sequence, quality, new_quality, trimmed_junk, junk_quality, has_adapter_file);

		/* write debugging files to outputs */
		if (read.length() >= minimum_read_length) {

			// add into hash table
			string key = get_key(info);
			DNAread tempLine;
			tempLine.sequence = read;
			tempLine.quality = new_quality;
			tempLine.coor = key;
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


string concat_reads (string sequence1, string sequence2, string quality1, string quality2, string &new_quality) {
	string concat = "";
	int highest_score = 0;

	// check one direction
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



void MPI_trim_and_match (char *infile, char *outfile, bool has_adapter_file) {

	ifstream in(infile);

	int size, rank, counter = 0;

	MPI_init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);


	// output filename
	stringstream ss;
	ss << outfile << "_" << rank << ".fastq" << endl;
	string out_filename(ss.str());
	ofstream out(out_filename.c_str());

	string info;
	while (getline(in, info)) {
		
		string sequence, na, quality;
		string read, new_quality, trimmed_junk, junk_quality;
		
		// read lines from in1
		getline(in, sequence);
		getline(in, na);
		getline(in, quality);

		// process only if line is assigned
		if (counter % size == rank) {

			//get key for hash table
			string key = get_key(info);

			//strip T's for file 1
			strip_t(sequence, quality);

			// trim read 1
			read = trim_read (sequence, quality, new_quality, trimmed_junk, junk_quality, has_adapter_file);

			if (read.length() >= minimum_read_length) {
				// check if key exists
				unordered_map<string, DNAread>::const_iterator found = myMap.find(key);

				//if no such key found
				if (found == myMap.end()) {

					// check quality
					if (quality_check(new_quality)) {
						out << "read2:" << key << endl
							<< read << endl
							<< na << endl
							<< new_quality << endl;
					}

				} else {
					// concatenate
					string concat_read, concat_quality;
					string read1 = myMap[key].sequence;
					string quality1 = myMap[key].quality;


					// match hash table key
					concat_read = concat_reads(read1, read, quality1, new_quality, concat_quality);
					
					if (concat_read.length() > 0) {

						// check quality
						if (quality_check(concat_quality)) {
							out << "concatenated:" << key << endl
								<< concat_read << endl
								<< na << endl
								<< concat_quality << endl;
						}

					} else {

						// check quality1
						if (quality_check(quality1)) {
							out << "read1:" << key << endl
								<< read1 << endl
								<< na << endl
								<< quality1 << endl;
						}

						// check quality2
						if (quality_check(new_quality)) {
							out << "read2:" << key << endl
								<< read << endl
								<< na << endl
								<< new_quality << endl;
						}

					}
				}
			}

		} else {
			counter++;
		}
	}
	out.close();

	// concatenate files together
	if (rank == 0) concat_files(outfile, size);

	MPI_Finalize();
}

void concat_files (char *outfile, int size) {
	stringstream ss;
	ss << outfile << ".fastq";
	string f = ss.str();
	ofstream out(f.c_str());
	for (int i = 0; i < size; i++) {
		stringstream t;
		t << outfile << "_" << i << ".fastq";
		string temp = t.str();
		ifstream in(temp.c_str());
		out << in.rdbuf();
		in.close();
	}
	out.close();
}


// single threaded read processing function
void process_reads (char* infile1, char* infile2, char* outfile, bool has_adapter_file) {

	int discarded1 = 0, discarded2 = 0, num_concat = 0, num_final = 0;

	// files for infile1
	ifstream in1(infile1);
	string out_filename1 = "trimmed_1_" + string(outfile);
	ofstream out1(out_filename1.c_str());
	string junk_filename1 = "junk_1_" + string(outfile);
	ofstream junk_out1(junk_filename1.c_str());
	string discarded_filename1 = "discarded_1_" + string(outfile);
	ofstream discarded_out1(discarded_filename1.c_str());
	string junk_discarded_filename1 = "junk_discarded_1_" + string(outfile);
	ofstream junk_discarded_out1(junk_discarded_filename1.c_str());


	// files for infile2
	ifstream in2(infile2);
	string out_filename2 = "trimmed_2_" + string(outfile);
	ofstream out2(out_filename2.c_str());
	string junk_filename2 = "junk_2_" + string(outfile);
	ofstream junk_out2(junk_filename2.c_str());
	string discarded_filename2 = "discarded_2_" + string(outfile);
	ofstream discarded_out2(discarded_filename2.c_str());
	string junk_discarded_filename2 = "junk_discarded_2_" + string(outfile);
	ofstream junk_discarded_out2(junk_discarded_filename2.c_str());


	// file for concatenated reads
	string concat_out_filename = "concat_" + string(outfile);
	ofstream concat_out(concat_out_filename.c_str());
	string final_out_filename = "final_" + string(outfile);
	ofstream final_out(final_out_filename.c_str());


	// trim fastq file
	// read 4 lines in a roll (4 lines = 1 set of data)
	// for (int i = start; i < end; i++) {

	//start clock
	clock_t start_time = clock();

	int j = 0;
	string info1, info2;
	while (getline(in1, info1) && getline(in2, info2)) {

		if (j % 100 == 0) {
			cout << "Reading line " << j+1 << "..." << endl;
		}
		
		j++;
		
		string sequence1, extra1, quality1, sequence2, extra2, quality2;
		string read1, trimmed_junk1, junk_quality1, new_quality1, read2, trimmed_junk2, junk_quality2, new_quality2;

		// read lines from the input file
		// file 1
		// getline(in, info); // get fastq read info
		getline(in1, sequence1); // get read 1
		getline(in1, extra1); // N/A
		getline(in1, quality1); // get the quality scores 1

		// file 2
		getline(in2, sequence2); // get read 2
		getline(in2, extra2); // N/A
		getline(in2, quality2); // get the quality score 2


		//strip T's for file 1
		strip_t(sequence1, quality1);

		//strip T's for file 1
		strip_t(sequence2, quality2);
		
		// trim read 1
		read1 = trim_read (sequence1, quality1, new_quality1, trimmed_junk1, junk_quality1, has_adapter_file);

		// trim read 2
		read2 = trim_read (sequence2, quality2, new_quality2, trimmed_junk2, junk_quality2, has_adapter_file);


		/* write debugging files to outputs */
		if (read1.length() >= minimum_read_length) {
			// write to output file
			out1 << info1 << endl
				<< read1 << endl
				<< extra1 << endl
				<< new_quality1 << endl;

			// write to junk output file
			junk_out1 << info1 << endl
					 << trimmed_junk1 << endl
					 << extra1 << endl
					 << junk_quality1 << endl;
		} else {
			discarded1++;
			// write to output file
			discarded_out1 << info1 << endl
						  << read1 << endl
						  << extra1 << endl
						  << new_quality1 << endl;

			// write to junk output file
			junk_discarded_out1 << info1 << endl
							   << trimmed_junk1 << endl
							   << extra1 << endl
							   << junk_quality1 << endl;
		}


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

		bool has_concat = false;
		string concatenated, new_quality;

		// concatenate reads
		if (read1.length() >= minimum_read_length && read2.length() >= minimum_read_length) {

			Coordinate coord1, coord2;
			coord1 = find_coord(info1);
			coord2 = find_coord(info2);

			if (is_same_coord(coord1,coord2)) {
				// if ((num_concat) % 100 == 0) {
				// 	cout << "read 1: (" << coord1.x << "," << coord1.y << "), "
				// 		 << "read 2: (" << coord2.x << "," << coord2.y << ")" << endl;
				// }

				concatenated = concat_reads(read1, read2, quality1, quality2, new_quality);

				if (concatenated.length() > 0) {
					num_concat++;
					has_concat = true;

					concat_out << info1 << endl
							   << concatenated << endl
							   << extra1 << endl
							   << new_quality << endl;
				}
			} else {
				cout << "Different coordinates; concatenation skipped {"
					 << "read 1: (" << coord1.x << "," << coord1.y << "), "
					 << "read 2: (" << coord2.x << "," << coord2.y << ")" 
					 << "}...";
			}
		}

		if (has_concat) {
			bool pass_quality = quality_check(new_quality);
			if (pass_quality) {
				num_final++;

				final_out << info1 << endl
						  << concatenated << endl
						  << extra1 << endl
						  << new_quality << endl;
			}
		}

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
	cout << "NUmber of reads that passed quality check: " << num_final << endl;


	in1.close();
	out1.close();
	junk_out1.close();
	discarded_out1.close();
	junk_discarded_out1.close();


	in2.close();
	out2.close();
	junk_out2.close();
	discarded_out2.close();
	junk_discarded_out2.close();

	concat_out.close();
	final_out.close();

	return;
}

void print_argument_directions(char* program) {
	cerr << "Usage: " << program << "<path_to_file1> <path_to_file2> <output_filename> " << endl;
	cerr << "optional: --adapter <path_to_adapter_file> (default: truseq adapters)";
	cerr << "             (*note that using your own adapter file may slow the program down)" << endl;
	cerr << "          --trim <minimum length of trimmed reads> (default: 11)" << endl;
	cerr << "          --match <minimum match length> (default: 7)" << endl;
}


int main (int argc, char** argv) {

	if (argc < 4) {
		print_argument_directions(argv[0]);
		exit(1); 
	}

	bool using_adapter_file = false;

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
			} else {
				print_argument_directions(argv[0]);
				exit(0);
			}

			i++;
		}
	}


	return 0;
}



