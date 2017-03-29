
#include <cstdlib>
#include <cstdio>
#include <string.h>
#include <string>
#include <vector>
#include <iomanip>
#include <iostream>
#include <fstream>
#include <mpi.h>
#include "global.hpp"
#include "MPI_readFastq.hpp"
#include "single_readFastq.hpp"
#include "memory_usage.hpp"
// #include "adapters.hpp"

using namespace std;


void print_argument_directions(char* program) {
	cerr << "Usage: " << program << " <path_to_file1> <path_to_file2> <output_filename> " << endl;
	cerr << "optional: --adapter <path_to_adapter_file> (default: truseq adapters)";
	cerr << "             (*note that using your own adapter file may slow the program down)" << endl;
	cerr << "          --trim <minimum length of trimmed reads> (default: 11)" << endl;
	cerr << "          --match <minimum match length> (default: 7)" << endl;
	cerr << "          --MPI (flag for using MPI)" << endl;
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

// guesstimate if the file is compressed by looking the file extension
int check_compression (string input) {
	vector<string> compression_exts = {"gz", "zip", "rar", "7z", "tar", "tgz", "tbz2", "bz2", "Z"};

	string ext = input;
	size_t pos;
	while (pos = ext.find(".") != string::npos) {	// look for the last "."
		ext = ext.substr(pos);
	}

	for (int i = 0; i < compression_exts.size(); i++) {
		if (ext.compare(compression_exts[i]) == 0) {
			return 1;
		}
	}

	return 0;
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
			} else if (string(argv[i]).compare("--MPI") == 0) {
				use_MPI = true;
			} else if (string(argv[i]).compare("--debug") == 0) {
				debug = true;
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

		MPI_process_reads(argv[1], argv[2], argv[3], rank, size, adapters, debug);

		MPI_Finalize();
	} else {
		process_reads(argv[1], argv[2], argv[3], 2, adapters, debug);
	}

	return 0;
}
