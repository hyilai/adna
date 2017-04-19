#ifndef GLOBAL
#define GLOBAL

// defaults
#define DEFAULT_MIN_CUT_LENGTH 11
#define DEFAULT_MIN_MATCH_LENGTH 7
#define DATABLOCKS 10000000
#define LINEBLOCKS 1000

// // universal sequence for tru seq adapters
// extern const std::string TRU_SEQ_ADAPTER = "GATCGGAAGAGCACACGTCTGAACTCCAGTCAC";

// // adapters
// extern std::vector<std::string> adapters;


// struct for holding info for a read
struct DNAread {
	// std::string info;
	std::string sequence;
	std::string quality;
};

// info 
// example: @HWI-D00550:263:C6V0DANXX:8:1101:1244:2115 1:N:0:CAGATC
struct DNAinfo {
	std::string first;	// first part of the info without the coordinate
	std::string second;	// second part of the info
};


// cutoff values
extern int minimum_read_length;
extern int minimum_match_length;

// set functions
void set_read_length(int);
void set_match_length(int);
// void set_adapters(char*);

#endif