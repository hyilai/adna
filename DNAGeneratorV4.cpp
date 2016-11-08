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

	if(argc != 5) {
		cout << "Error: Must provide arguments <length of DNA> <length of adapters> <number of adapters> <number of fragments>" << endl;
		return 0;
	}

	int Length_DNA = atoi(argv[1]); //This gets the length of the DNA from the command argument
	int Length_adapter = atoi(argv[2]); //This gets the length of the adapter from the command argument 
	int DNA_split = atoi(argv[3]); //This is how many times we will split the DNA  
	int number_of_adapters = atoi(argv[4]); //This is the amount of adapters we have
	string DNA;
	string adapter;
	string adapters[number_of_adapters];

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


	for(int k = 0; k < number_of_adapters; k++) {
		//This randomly generates an adapter
		for(i=0;i<Length_adapter;i++) {
			ss << get_random_base();
		}

		ss >> adapter;

		adapters[k] = adapter;

		ss.str("");
		ss.clear();
		adapter = "";

	}


	string DNA_adapters;
	ofstream outFile("Test.fastq");
	ofstream out("testdna.txt");
	ofstream extraDNAFile("extra.txt");

	int len = Length_DNA/DNA_split;

	//Split DNA and put adapters in between
	for(i=0;i<DNA_split;i++) {

		outFile << endl;
		int start = i*(Length_DNA/DNA_split);
		int end = (i+1)*(Length_DNA/DNA_split);



		//Generates extra DNA for the ends of the adapters
		//The length of the extra DNA will be between 0 - 3/4(Length_adapter) 
		string extra1 = make_extra_DNA((Length_adapter*3)/4);
		string extra2 = make_extra_DNA((Length_adapter*3)/4);



		//This randomly selects an adapter from the adapters array
		int index =  rand() % number_of_adapters;
		adapter = adapters[index];

		// cout << adapter << " " << index << " " << number_of_adapters << endl;




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


		extraDNAFile << extra1 << "  " << adapterFrag1 << "  " << DNA.substr(start, len);
		if (i == DNA_split-1 && end != Length_DNA) {
			extraDNAFile << DNA.substr(end, Length_DNA);
		}
		extraDNAFile << "  " << adapterFrag2 << "  " << extra2 << endl;
		

		outFile << adapterFrag2 << extra2 << endl;
		outFile << "+" << endl;
		

		int line_length = extra1.length() + extra2.length() + adapterFrag1.length() + adapterFrag2.length() + len;

		for(int k=0;k<line_length;k++) {
			outFile << "Q" ;
		}

		outFile <<  endl;


		out << DNA.substr(start, len) << endl;
	}

	extraDNAFile.close();
	outFile.close();
	out.close();

	ofstream adapterFile("adapter.fa");

	for(int k = 0; k<number_of_adapters;k++) {
		adapterFile << ">TestAdapter" << k+1 << endl;
		adapterFile << adapters[k] << endl;
	}

	

	adapterFile.close();

	//Prints results
	//cout << "DNA:" << DNA << endl;
	//cout << "Adapter:"<< adapter << endl;
}
