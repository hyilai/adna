#ifndef SINGLE_READ
#define SINGLE_READ

struct DNAread;

// Hash table
// std::unordered_map<std::string, DNAread> myMap;


void trim_file(char*, char*, int, bool, int&);
void process_reads (char*, char*, char*, int, bool);


#endif