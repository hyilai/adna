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
		cout << "Error: Must provide arguments <length of DNA> <number of fragments>" << endl;
		return 0;
	}

	int Length_DNA = atoi(argv[1]); //This gets the length of the DNA from the command argument
	int DNA_split = atoi(argv[2]); //This is how many times we will split the DNA 
	string DNA;
	vector<string> adapters;

	// get fasta file
	ifstream a_file("adapters.fasta");
	string line;
	while(getline(a_file, line)) {
		getline(a_file, line);
		adapters.push_back(line);
	}

	int number_of_adapters = adapters.size();

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


	string DNA_adapters;
	ofstream outFile("Test.fastq");
	ofstream out("testdna.txt");
	ofstream extraDNAFile("extra.txt");

	int len = Length_DNA/DNA_split;

	//Split DNA and put adapters in between
	for(i=0;i<DNA_split;i++) {

		outFile << i+1 << endl;
		int start = i*(Length_DNA/DNA_split);
		int end = (i+1)*(Length_DNA/DNA_split);


		//This randomly selects an adapter from the adapters array
		int index =  rand() % number_of_adapters;
		string adapter = adapters[index];


		//Generates extra DNA for the ends of the adapters
		//The length of the extra DNA will be between 0 - 3/4(Length_adapter)
	
		int Length_adapter = adapter.length();

		string extra1 = make_extra_DNA((Length_adapter*3)/4);
		string extra2 = make_extra_DNA((Length_adapter*3)/4);



		//cout << start << " " << end << "  " << DNA.substr(start, len) << endl;

		// cut a chunk of the adapter out
		// length of adapter fragment is between 7 to length of adapter
		// randomly add bases to the end of adapter
		int temp = rand() % (Length_adapter - MIN_ADAPTER_LENGTH) + MIN_ADAPTER_LENGTH;
		// cout << adapter.length() << " " << Length_adapter << " " << Length_adapter-temp << endl;
		string adapterFrag1 = adapter.substr(Length_adapter-temp,temp);
		string adapterFrag2 = adapter.substr(0,temp);


		int num_extra = rand() % 15;

		for (int k = 0; k < num_extra; k++) {
			int b = rand() % 4;

		}

		//outFile << extra1 << adapterFrag1;
		outFile << DNA.substr(start, len);
		if (i == DNA_split-1 && end != Length_DNA) {
			outFile << DNA.substr(end, Length_DNA);
		}


		//extraDNAFile << extra1 << "  " << adapterFrag1 << "  " << DNA.substr(start, len);
		extraDNAFile << DNA.substr(start, len);
		if (i == DNA_split-1 && end != Length_DNA) {
			extraDNAFile << DNA.substr(end, Length_DNA);
		}
		extraDNAFile << "  " << adapterFrag2 << "  " << extra2 << endl;
		

		outFile << adapterFrag2 << extra2 << endl;
		outFile << "+" << endl;
		

		//int line_length = extra1.length() + extra2.length() + adapterFrag1.length() + adapterFrag2.length() + len;
		int line_length =  extra2.length() + adapterFrag2.length() + len;


		for(int k=0;k<line_length;k++) {
			outFile << "Q" ;
		}

		outFile <<  endl;


		out << DNA.substr(start, len) << endl;
	}

	extraDNAFile.close();
	outFile.close();
	out.close();

	//Prints results
	//cout << "DNA:" << DNA << endl;
	//cout << "Adapter:"<< adapter << endl;
}
