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
#include <mpi.h>
#include "SmithWaterman.hpp"

using namespace std;

class SmithWaterman;

// containers
vector<vector<string> > file_line;
vector< vector<string> > trimmed_seq;
vector<string> adapters;

void parse_file (int start, int end) {

	// trim fastq file
	// read 4 lines in a roll (4 lines = 1 set of data)
	for (int i = start; i < end; i++) {

		cout << "Parsing line " << i << "..." << endl;
		string line1, line2, line3, line4, trimmed, quality;

		line1 = file_line[i][0];
		line2 = file_line[i][1];
		line3 = file_line[i][2];
		line4 = file_line[i][3];

		for (int i = 0; i < adapters.size(); i++) {

			SmithWaterman* sw = new SmithWaterman(line2,adapters[i], line4, "", 0);

			// string temp = sw->trim_both_sides();	// trim read

			if (temp.length() > 0 && temp.length() < line2.length()) {
				trimmed = temp;
				quality = sw->get_quality1();	// get trimmed quality control string
			}

			trimmed_seq[j].push_back(temp);

			delete sw;
		}

		// output new data
		//out << line1 << endl
		//	<< trimmed << endl
		//	<< line3 << endl
		//	<< quality << endl;
		// WIP; put them into the vector for now
		file_line[i][1] = trimmed;
		file_line[i][3] = quality;
	}

	file.close();
	out.close();

	return 0;
}


int main (int argc, char** argv) {

	if (argc != 4) {
		cerr << "Arguments: <fastq_file> <adapter_file> <output_file>" << endl;
		exit(1); 
	}

	ifstream file(argv[1]);
	ifstream adpt(argv[2]);

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

	// get all file lines
	string line;
	while (getline(file,line)) {
		vector<string> line_group;
		line_group.push_back(line);	// first line; information of the read
		getline(file, line);
		line_group.push_back(line);	// second line; sequence read
		getline(file, line);
		line_group.push_back(line);	// third line
		getline(file, line);
		line_group.push_back(line);	// fourth line; quality control
		file_lines.push_back(line_group);
	}

	ofstream out(argv[3]);

	// set up MPI stuff
	int comm_sz;
	int my_rank;

	MPI_Init(NULL,NULL);
	MPI_comm_size(MPI_COMM_WORLD, &comm_sz);
	MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

	if (my_rank != 0) {
		//MPI_Send();
	} else {
		//MPI_Recv();
	}

	MPI_Finalize();

	return 0;
}



