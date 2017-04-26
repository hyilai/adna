
#include <cstdlib>
#include <cstdio>
#include <string.h>
#include <string>
#include <sstream>
#include <vector>
#include <iomanip>
#include <iostream>
#include <fstream>
#include <mpi.h>
#include "global.hpp"
#include "MPI_readFastq.hpp"
// #include "single_readFastq.hpp"
#include "memory_usage.hpp"

using namespace std;


void print_argument_directions(char* program) {
	cerr << "Usage: " << program << " <path_to_file1> <path_to_file2> <output_filename> " << endl;
	cerr << "optional: --adapter <path_to_adapter_file> (default: truseq adapters)";
	cerr << "             (*note that using your own adapter file may slow the program down)" << endl;
	cerr << "          --trim <minimum length of trimmed reads> (default: 11)" << endl;
	cerr << "          --match <minimum match length> (default: 7)" << endl;
	// cerr << "          --MPI (flag for using MPI)" << endl;
}


vector<string> set_adapters (char* c_filename) {
	// get the list of adapters
	ifstream adpt(c_filename);
	string ad;
	vector<string> adapters;
	while (1) {
		// if (adpt.peek() == EOF) {
		// 	break;
		// }

		//discard first line
		if (!getline(adpt, ad)) break;

		//get adapter sequence in the second line
		getline(adpt, ad);
		adapters.push_back(ad);
	}
	adpt.close();
	return adapters;
}

// guesstimate if the file is a gzipped file by looking the file extension
bool check_compression (string input) {
	// vector<string> compression_exts = {"gz", "zip", "rar", "7z", "tar", "tgz", "tbz2", "bz2", "Z"};

	// string ext = input;
	// size_t pos;
	// while (pos = ext.find(".") != string::npos) {	// look for the last "."
	// 	ext = ext.substr(pos);
	// }

	// for (int i = 0; i < compression_exts.size(); i++) {
	// 	if (ext.compare(compression_exts[i]) == 0) {
	// 		return 1;
	// 	}
	// }

	if (input.find(".gz") != string::npos) {
		return true;
	}

	return false;
}

string get_filename(char* filepath) {
	// erase path
	string temp(filepath);
	size_t found = temp.find_last_of("/");
	if (found != string::npos) {
		temp = temp.substr(found+1);
	}
	// erase .gz extension
	found = temp.find(".gz");
	if (found != string::npos) {
		temp = temp.substr(0,found);
	}
	return temp.c_str();
}


// check if the two input files are correctly paired
bool match_files (string &file1, string &file2) {
	string line1, line2;

	// get the first line in file1
	ifstream in1(file1.c_str());
	getline(in1, line1);
	in1.close();

	// get the first line in file2
	ifstream in2(file2.c_str());
	getline(in2, line2);
	in2.close();

	string line1_first, line1_second, line2_first, line2_second;
	int line1_num, line2_num;

	// grab the info from file1
	size_t found = line1.find(" ");
	if (found) {
		line1_first = line1.substr(0, found);
		line1_second = line1.substr(found+1);

		found = line1_first.find_last_of(":");			// y coordinate
		line1_first = line1_first.substr(0, found);
		found = line1_first.find_last_of(":");			// x coordinate
		line1_first = line1_first.substr(0, found);
		found = line1_first.find_last_of(":");			// flowcell
		line1_first = line1_first.substr(0, found);
		
		found = line1_second.find(":");					// read number (1 or 2)
		line1_num = atoi(line1_second.substr(0, found).c_str());
		line1_second = line1_second.substr(found+1);
	}

	// grab the info from file2
	found = line2.find(" ");
	if (found) {
		line2_first = line2.substr(0, found);
		line2_second = line2.substr(found+1);

		found = line2_first.find_last_of(":");			// y coordinate
		line2_first = line2_first.substr(0, found);
		found = line2_first.find_last_of(":");			// x coordinate
		line2_first = line2_first.substr(0, found);
		found = line2_first.find_last_of(":");			// flowcell
		line2_first = line2_first.substr(0, found);
		
		found = line2_second.find(":");					// read number (1 or 2)
		line2_num = atoi(line2_second.substr(0, found).c_str());
		line2_second = line2_second.substr(found+1);
	}

	// cout << line1_num << " " << line2_num << endl;

	// check if the files are of the same pair reads (based on the info grabbed from the files)
	if (line2_first.compare(line1_first) == 0) {

		// switch the filenames around if the filenames are not given in order
		if (line1_num == 2 && line2_num == 1) {
			string temp = file1;
			file1 = file2;
			file2 = temp;
		}

		return true;
	}
	return false;
}


int main (int argc, char** argv) {


	if (argc < 4) {
		print_argument_directions(argv[0]);
		exit(1); 
	}

	//Initialize data and files to track RAM and CPU usage
	initCPUCurrent();
	initCPUCurrentProcess();

	// initialize values
	set_read_length(DEFAULT_MIN_CUT_LENGTH);
	set_match_length(DEFAULT_MIN_MATCH_LENGTH);


	bool using_adapter_file = false;
	bool use_MPI = false;
	bool debug = false;

	// universal sequence for tru seq adapters
	const string TRU_SEQ_ADAPTER = "GATCGGAAGAGCACACGTCTGAACTCCAGTCAC";

	// adapters
	vector<std::string> adapters;
	adapters.push_back(TRU_SEQ_ADAPTER);


	// grab values from arguments
	if (argc > 4) {
		int i = 4;
		while (i < argc) {

			if (string(argv[i]).compare("--adapter") == 0) {
				using_adapter_file = true;
				// add adapters
				adapters = set_adapters(argv[++i]);
			} else if (string(argv[i]).compare("--trim") == 0) {
				set_read_length(atoi(argv[++i]));
			} else if (string(argv[i]).compare("--match") == 0) {
				set_match_length(atoi(argv[++i]));
			} else if (string(argv[i]).compare("--debug") == 0) {
				debug = true;
			} else {
				print_argument_directions(argv[0]);
				exit(0);
			}

			i++;
		}
	}


	string file1(argv[1]);
	string file2(argv[2]);
	bool is_file1_gz = check_compression(file1);
	bool is_file2_gz = check_compression(file2);


	// unzip gz files to temp files
	if (is_file1_gz) {
		// get filenames from filepath
		file1 = get_filename(argv[1]);
		stringstream command;
		command << "gunzip -c " << argv[1] << " > " << file1;
		system(command.str().c_str());
	}

	if (is_file2_gz) {
		file2 = get_filename(argv[2]);
		stringstream command;
		command << "gunzip -c " << argv[2] << " > " << file2;
		system(command.str().c_str());
	}


	// initialize MPI
	int size, rank;
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);


	bool paired_file = match_files(file1, file2);

	// change const char* to char*
	char *infile1 = new char[file1.size() + 1];
	strcpy(infile1, file1.c_str());
	char *infile2 = new char[file2.size() + 1];
	strcpy(infile2, file2.c_str());


	// wait til all processes are done up to this point
	MPI_Barrier(MPI_COMM_WORLD);

	if (paired_file) {

		process_reads(infile1, infile2, argv[3], rank, size, adapters, debug);

	} else {
		if (rank == 0) cerr << "Error: input files are not paired" << endl;
	}

	if (rank == 0) {
		// remove temporary files
		if (is_file1_gz) remove(infile1);
		if (is_file2_gz) remove(infile2);
	}

	delete[] infile1;
	delete[] infile2;

	MPI_Finalize();


	return 0;
}
