#include <cstdlib>
#include <cstdio>
#include <string.h>
#include <vector>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <fstream>
#include <ctime>
#include "SmithWaterman.hpp"
#include <mpi.h>
#include <unordered_map>x
// #include <pthread.h>

using namespace std;

#define MIN_ADAPTER_LENGTH 7


class SmithWaterman;

struct Coordinate {
	int x, y;
};

struct DNAread {
	string sequence;
	string quality;
	string coor;
};


//Hash table 
unordered_map<string, DNAread> myMap;


int main (int argc, char** argv) {

	srand(time(NULL));
	
	int comm_sz;    /* number of processes  */
	int my_rank;    /* process rank         */
	
	
	if (argc != 3) {
		cerr << "Arguments: <minimum matching score> <output file>" << endl;
		exit(1); 
	}

	

	ifstream in1("testRead1.fastq");
	string info1, sequence1, na1, quality1;
	while(getline(in1, info1)) {

		// read lines from in1
		stringstream ss;
		ss.str(info1);
		string token;
		for (int i = 0; i < 2; i++) {
			getline(ss, token, ':');
		}
		
		
		string x;
		getline(ss, x, ':');
		string y;
		getline(ss, y, ':');
		
		string key = x + ":" + y;


		getline(in1, sequence1);
		getline(in1, na1);
		getline(in1, quality1);
		
		DNAread tempLine;
		
		tempLine.sequence = sequence1;
		tempLine.quality = quality1;
		tempLine.coor = key;
		
		myMap.insert({key,tempLine});
	}

	
	
	int match_score = atoi(argv[1]);

	ofstream out(argv[2]);	// result file
	
	ifstream in2("testRead2.fastq");

	MPI_Status status;
	MPI_File fh;
	MPI_offset offset;



	MPI_File_open(MPI_COMM_WORLD, argv[2], MPI_MODE_CREATE|MPI_MODE_WRONLY, MPI_INFO_NULL, &fh); 
	nints


	// containers
	vector<string> matched_seq;

	// analyze sw
	int j = 0;
	vector<double> total_elapsed_time;

	// clock start
	clock_t start = clock();
	
	//

	//Initialize the command line arguments on every process
	MPI_Init(&argc, &argv);
	
	//Get the number of processes in MPI_COMM_WORLD, and put 
	//it in the 'comm_sz" variable
	MPI_Comm_size(MPI_COMM_WORLD, &comm_sz); 
	MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);


	int counter = 0;
	
	string info2;
	while(getline(in2, info2)) {
		
		

		if((counter % comm_sz) ==  my_rank) {

			//
			string sequence2, na2, quality2;
			
			// read lines from in1
			getline(in2, sequence2);
			getline(in2, na2);
			getline(in2, quality2);

			// cout << "Reading line " << j+1 << "..." << endl;
			
			// read lines from in1
			stringstream ss;
			ss.str(info2);
			string token;
			for (int i = 0; i < 2; i++) {
				getline(ss, token, ':');
			}


			string x;
			getline(ss, x, ':');
			string y;
			getline(ss, y, ':');

			string key = x + ":" + y;

			
			string sequence1 = myMap[key].sequence;
			string quality1 = myMap[key].quality;
		

			if(key.compare(myMap[key].coor)==0) {
				string matched = "";
				int highest_score = 0;

				// start clock for sw algorithm
				clock_t start_time = clock();

				// check one direction
				SmithWaterman* sw = new SmithWaterman(sequence1, sequence2, quality1, quality2, match_score);
				int curr_high_score = sw->get_highest_score();
				if (highest_score < curr_high_score) {
					if (sw->match_reads()) {
						matched = sw->get_matched_string();
					}
				}
				delete sw;

				// end clock for sw algorithm
				clock_t end_time = clock();

				// for calculations
				total_elapsed_time.push_back(double(end_time - start_time) / CLOCKS_PER_SEC);

				// output file
				if (matched.length() > 0) {
					out << matched << endl;
					cout << matched << endl;
					MPI_File_write(fh, buf, 1000, MPI_INT, MPI_STATUS_IGNORE); 
				}
				
// 				DNAread* del = myMap[key];
// 				delete rel;
// 				myMap.erase(key);
			}
				
			
		}
		else {

			string sequence2, na2, quality2;

			// read lines from in2
			getline(in2, sequence2);
			getline(in2, na2);
			getline(in2, quality2);
		}
		j++;
		counter++;
	}
	
	MPI_File_close(&fh); 
	/*Inform MPI that this process has finished*/
	MPI_Finalize();
	
	
	// clock ends
	clock_t end = clock();
	double total_time_passed = double(end - start) / CLOCKS_PER_SEC;

	out.close();


	return 0;
}
