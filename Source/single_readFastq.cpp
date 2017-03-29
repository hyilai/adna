
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
#include "memory_usage.hpp"
#include "steps.hpp"
#include "utilities.hpp"
// #include "hash_table.hpp"
#include "global.hpp"
#include "single_readFastq.hpp"
#include "ZipLib/ZipFile.h"
#include "ZipLib/streams/memstream.h"
#include "ZipLib/methods/Bzip2Method.h"


using namespace std;

// Hash table
unordered_map<string, DNAread> myMap;


// trim one read file and put the data into the hash table; single-threaded
void trim_file (char *infile, char *outfile, int file_num, bool has_adapter_file, int &discarded) {

	discarded = 0;

	// ifstream in(infile);
	// decompression
	ZipArchive::Ptr archive = ZipFile::Open(string(infile));
	ZipArchiveEntry::Ptr entry = archive->GetEntry(0);
	istream* in = entry->GetDecompressionStream();

	// string line;
	// getline(**decompressStream, info);

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
	while (getline(*in, info)) {

		if (j++ % LINEBLOCKS == 0) {
			cout << "\tParsing read " << j << "..." << endl;
		}

		string sequence, extra, quality;
		string read, trimmed_junk, junk_quality, new_quality;

		// read lines from the input file
		// getline(*in, info);		// get fastq read info
		getline(*in, sequence);	// get read
		getline(*in, extra);		// get the third line
		getline(*in, quality);	// get the quality scores 


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



// single threaded read processing function
void process_reads (char* infile1, char* infile2, char* outfile, int file_num, bool has_adapter_file) {

	int discarded1 = 0, discarded2 = 0, num_concat = 0, num_final = 0;

	//start clock
	clock_t start_time = clock();

	cout << "Reading file1..." << endl;

	// get hash table for read file1
	trim_file(infile1, outfile, 1, has_adapter_file, discarded1);

	int num_reads1 = myMap.size();

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

	// ifstream in2(infile2);
	// decompression
	ZipArchive::Ptr archive = ZipFile::Open(string(infile)2);
	ZipArchiveEntry::Ptr entry = archive->GetEntry(0);
	istream* decompressStream = entry->GetDecompressionStream();

	int j = 0;
	string info2;
	while (getline(*in2, info2)) {

		if (j++ % LINEBLOCKS == 0) {
			cout << "\tParsing read " << j << "..." << endl;
		}
		
		string sequence2, extra2, quality2;
		string read2, trimmed_junk2, junk_quality2, new_quality2;

		// file 2
		getline(*in2, sequence2); // get read 2
		getline(*in2, extra2); // N/A
		getline(*in2, quality2); // get the quality score 2

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

			// get the key for read2
			string key = get_key(info2);

			// check if the hash table has said key
			if (myMap.find(key) == myMap.end()) {

				if (quality_check(read2)) {

					// output read2 on its own
					final_out << info2 << endl
							  << read2 << endl
							  << extra2 << endl
							  << new_quality2 << endl;
				}

			} else {

				string info1, read1, quality1, concatenated, concat_quality;

				info1 = myMap[key].info;
				read1 = myMap[key].sequence;
				quality1 = myMap[key].quality;

				myMap.erase(key);

				// concatenate reads
				concatenated = concat_reads(read1, read2, quality1, quality2, concat_quality);

				if (concatenated.length() > 0) {
					num_concat++;

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
					if (quality_check(quality1)) {
						final_out << info1 << endl
								  << read1 << endl
								  << extra2 << endl
								  << quality1 << endl;
					}

					if (quality_check(quality2)) {
						final_out << info2 << endl
								  << read2 << endl
								  << extra2 << endl
								  << new_quality2 << endl;
					}
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

	print_diagnostics(elapsed_time, num_reads1, j, discarded1, discarded2, num_concat, num_final);

	in2.close();
	out2.close();
	junk_out2.close();
	discarded_out2.close();
	junk_discarded_out2.close();

	final_out.close();
}
