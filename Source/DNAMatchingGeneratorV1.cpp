#include <iostream>
#include <stdio.h>
#include <vector>
#include <stdlib.h>
#include <cstring>
#include <sstream>
#include <cstdlib>
#include <ctime>
#include <fstream>

using namespace std;

#define MIN_ADAPTER_LENGTH 7

char get_random_base() {
	int random = rand() % 20 + 1;
	if(random < 6) {
		return 'A';
	} 
	else if(random < 11) {
		return 'G';
	} 
	else if(random < 16) {
		return 'T';
	}
	else {
		return 'C';
	}
}


//Generates DNA of random length between 0 and range
string make_extra_DNA(int range) {
	int length = rand() % range + 1;

	string DNA;

	int i;
	stringstream ss;
	for(i=0;i<length;i++) {
		ss << get_random_base();
	}

	ss >> DNA;

	ss.str("");
	ss.clear();

	return DNA;

}

string make_DNA (int range) {
	string DNA;

	int i;
	stringstream ss;
	for(i=0;i<range;i++) {
		ss << get_random_base();
	}

	ss >> DNA;

	ss.str("");
	ss.clear();

	return DNA;

}


int main(int argc, char** argv) {

	if(argc != 3) {
		cout << "Error: Must provide arguments <number of fragments> <length of fragment>" << endl;
		return 0;
	}

	int fragments = atoi(argv[1]); //This gets the length of the DNA from the command argument
	int DNA_split = atoi(argv[2]); //This is how many times we will split the DNA 
	// string DNA;

	//This is for the random number generator
	srand(time(NULL));


	//This randomly generates the DNA
	// int i;
	// stringstream ss;
	// for(i=0;i<Length_DNA;i++) {
	// 	ss << get_random_base();
	// }

	// ss >> DNA;

	// ss.str("");
	// ss.clear();

	ofstream outFile1("testRead1.fastq");
	ofstream outFile2("testRead2.fastq");

	int i = 0;
	// outFile << DNA << endl;

	//Split DNA and put adapters in between
	while(i++ < fragments) {
		int overlap = rand() % (DNA_split/2);
		overlap = (overlap >= 7) ? overlap : 7;
		int temp1 = rand() % (DNA_split - overlap) + 1;
		int temp2 = DNA_split - temp1 - overlap;
		while (temp2 == 0) {
			temp1 = rand() % (DNA_split - overlap) + 1;
			temp2 = DNA_split - temp1 - overlap;
		}

		int actual_length = temp1 + temp2 + overlap;

		string DNA = make_DNA(actual_length);

		string read1 = DNA.substr(0, temp1 + overlap);

		string read2 = DNA.substr(temp1, temp2 + overlap);


		outFile1 << i << endl;
		outFile1 << read1 << endl;
		outFile1 << "" << endl;
		int j;
		for(j = 0; j<DNA_split;j++) {
			outFile1 << "+";
		} 
		outFile1 << endl;

		outFile2 << i << endl;
		outFile2 << read2 << endl;
		outFile2 << "" << endl;
		for(j = 0; j<DNA_split;j++) {
			outFile2 << "+";
		} 
		outFile2 << endl;

	}

	outFile1.close();
	outFile2.close();

	//Prints results
	//cout << "DNA:" << DNA << endl;
	//cout << "Adapter:"<< adapter << endl;
}
