#include <iostream>
#include <stdio.h>
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


int main(int argc, char** argv) {

	int Length_DNA = atoi(argv[1]); //This gets the length of the DNA from the command argument
	int Length_adapter = atoi(argv[2]); //This gets the length of the adapter from the command argument 
	int DNA_split = atoi(argv[3]); //This is how many times we will split the DNA  

	string DNA;
	string adapter;

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
	//This randomly generates the adapter
	for(i=0;i<Length_adapter;i++) {
		ss << get_random_base();
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

		// cut a chunk of the adapter out
		// length of adapter fragment is between 7 to length of adapter
		// randomly add bases to the end of adapter
		int temp = rand() % (Length_adapter-1) + MIN_ADAPTER_LENGTH;
		string adapterFrag1 = adapter.substr((Length_adapter-1)-temp,temp);
		string adapterFrag2 = adapter.substr(0,temp);

		stringstream temp2;
		temp2 << adapterFrag1;

		int num_extra = rand() % 15;

		for (int k = 0; k < num_extra; k++) {
			int b = rand() % 4;

		}

		outFile << adapterFrag1;
		outFile << DNA.substr(start, len);
		if (i == DNA_split-1 && end != Length_DNA) {
			outFile << DNA.substr(end, Length_DNA);
		}


		outFile << adapterFrag2 << endl;
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