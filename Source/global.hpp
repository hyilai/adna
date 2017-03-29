#ifndef GLOBAL
#define GLOBAL

// defaults
#define DEFAULT_MIN_CUT_LENGTH 11
#define DEFAULT_MIN_MATCH_LENGTH 7
#define LINEBLOCKS 500

// // universal sequence for tru seq adapters
// extern const std::string TRU_SEQ_ADAPTER = "GATCGGAAGAGCACACGTCTGAACTCCAGTCAC";

// // adapters
// extern std::vector<std::string> adapters;


// struct for holding info for a read
struct DNAread {
	std::string info;
	std::string sequence;
	std::string quality;
};


// cutoff values
extern int minimum_read_length;
extern int minimum_match_length;

// set functions
void set_read_length(int);
void set_match_length(int);
// void set_adapters(char*);

#endif