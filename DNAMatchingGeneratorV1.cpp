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


int main(int argc, char** argv) {

	if(argc != 3) {
		cout << "Error: Must provide arguments <length of DNA> <length of fragments>" << endl;
		return 0;
	}

	int Length_DNA = atoi(argv[1]); //This gets the length of the DNA from the command argument
	int DNA_split = atoi(argv[2]); //This is how many times we will split the DNA 
	string DNA;

	//This is for the random number generator
	srand(time(NULL));


	//This randomly generates the DNA
	int i;
	stringstream ss;
	for(i=0;i<Length_DNA;i++) {
		ss << get_random_base();
	}

	ss >> DNA;

	ss.str("");
	ss.clear();

	ofstream outFile1("testRead1.fastq");
	ofstream outFile2("testRead2.fastq");

	int prev_Location1 = 5;
	int prev_Location2 = 5;
	int counter = 0;
	i = 0;
	// outFile << DNA << endl;

	//Split DNA and put adapters in between
	while(counter < Length_DNA) {
		i++;
		int temp = rand() % (5);
		temp = temp + prev_Location1;
		counter = temp + DNA_split;
		if(counter > Length_DNA) {
			break;
		}

		int start1 = temp;
		string read1 = DNA.substr(start1,DNA_split);
		prev_Location1 = start1+DNA_split+5;


		outFile1 << i << endl;
		outFile1 << read1 << endl;
		outFile1 << "" << endl;
		int j;
		for(j = 0; j<DNA_split;j++) {
			outFile1 << "+";
		} 
		outFile1 << endl;

		
		temp = rand() % (5);
		temp = temp + DNA_split/2 + prev_Location2;
		counter = temp +  + DNA_split;
		if(counter > Length_DNA) {
			break;
		}

		int start2 = temp+1;
		string read2 = DNA.substr(start2,DNA_split);
		prev_Location2 = start2+DNA_split+5;

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
