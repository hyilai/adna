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

	if(argc != 3) {
		cout << "Error: Must provide arguments <total length of string> <length of each line>" << endl;
		return 0;
	}
	int length = atoi(argv[1]); //This gets the length of the random string from the command argument
	int length_line = atoi(argv[2]); //This gets the length of the random string from the command argument

	//Output File
	ofstream outFile("random.fastq");


	//This generates the random string and formats it to look like a fastq file 
	int i;
	int j;
	for(j=0;j<length;j++) {
		outFile <<  endl;
		for(i=0;i<length_line;i++,j++) {
			outFile << get_random_base();
		}
		outFile <<  endl;
		outFile << "+" << endl;
		outFile << "QQQQQQQQQQQQQ" << endl;
	}

	outFile.close();

}