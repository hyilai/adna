#include <cstdlib>
#include <cstdio>
#include <string.h>
#include <fstream>
#include <iostream>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/copy.hpp>
#include <boost/iostreams/filter/gzip.hpp>

using namespace std;

int main () {
	// decompression
	// igzstream in("archive.gz");
	// open and read from gzipped file
	ifstream file("/local_storage/capstone/15-006-001_ATCACG_L004_R1_001.fastq.gz", ios_base::in | ios_base::binary);
	boost::iostreams::filtering_istream in;
	in.push(boost::iostreams::gzip_decompressor());
	in.push(file);

	string line;
	getline(in, line);
	cout << line << endl;


	return 0;
}