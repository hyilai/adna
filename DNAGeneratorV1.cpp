#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <cstring>
#include <sstream>
#include <cstdlib>
#include <ctime>
#include <fstream>

using namespace std;

int main(int argc, char** argv) {

	int Length_DNA = atoi(argv[1]); //This gets the length of the DNA from the command argument
	int Length_adapter = atoi(argv[2]); //This gets the length of the adapter from the command argument 
	int DNA_split = atoi(argv[3]); //This is how many times we will split the DNA  

	string DNA;
	string adapter;

	//This is for the random number generator
	srand(time(NULL));


	//This randomly generates the DNA
	int i, random;
	stringstream ss;
	for(i=0;i<Length_DNA;i++) {

		random = rand() % 20 + 1;
		if(random < 6) {
			ss << 'A';
		} 
		else if(random < 11) {
			ss << 'G';
		} 
		else if(random < 16) {
			ss << 'T';
		}
		else if(random < 21) {
			ss << 'C';
		}
	}

	ss >> DNA;

	ss.str("");
	ss.clear();
	//This randomly generates the adapter
	for(i=0;i<Length_adapter;i++) {
		random = rand() % 20 + 1;
		if(random < 6) {
			ss << 'A';
		} 
		else if(random < 11) {
			ss << 'G';
		} 
		else if(random < 16) {
			ss << 'T';
		}
		else if(random < 21) {
			ss << 'C';
		}
	}

	ss >> adapter;


	string DNA_adapters;
	ofstream outFile("Test.fastq");
	ofstream out("testdna.txt");

	int len = Length_DNA/DNA_split;

	//Split DNA and put adapters in between
	for(i=0;i<DNA_split;i++) {

		outFile << endl;
		int start = i*(Length_DNA/DNA_split);
		int end = (i+1)*(Length_DNA/DNA_split);

		cout << start << " " << end << "  " << DNA.substr(start, len) << endl;

		outFile << adapter << DNA.substr(start, len);
		if (i == DNA_split-1 && end != Length_DNA) {
			outFile << DNA.substr(end, Length_DNA);
		}
		outFile << adapter << endl;
		outFile << "+" << endl;
		outFile << "QQQQQQQQQQQQQ" << endl;
		out << DNA.substr(start, len) << endl;
	}

	outFile.close();
	out.close();

	ofstream adapterFile("adapter.fa");

	adapterFile << ">TestAdapter1" << endl;
	adapterFile << adapter << endl;

	adapterFile.close();

	//Prints results
	cout << "DNA:" << DNA << endl;
	cout << "Adapter:"<< adapter << endl;
}