#include <cstdlib>
#include <cstdio>
#include <string.h>
#include <fstream>
#include <iostream>
// #include "ZipLib/ZipFile.h"
// #include "ZipLib/streams/memstream.h"
// #include "ZipLib/methods/Bzip2Method.h"
#include "zlib-1.2.11/zlib.h"
#include "zlib-1.2.11/gzstream.h"

using namespace std;

int main () {
	// decompression
	// ZipArchive::Ptr archive = ZipFile::Open("archive.gz");
	// if (archive == nullptr) {
	// 	cerr << "cant open gz file" << endl;
	// 	exit(1);
	// } else {
	// 	cout << archive->GetEntriesCount() << endl;
	// }
	// ZipArchiveEntry::Ptr entry = archive->GetEntry(0);
	// istream* in = entry->GetDecompressionStream();

	igzstream in("archive.gz");
	string line;
	getline(in, line);
	cout << line << endl;

	return 0;
}