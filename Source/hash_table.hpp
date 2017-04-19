#ifndef H_HASH_TABLE
#define H_HASH_TABLE

#include "global.hpp"
struct DNAread;

class hash_table {
private:
	// Hash table

public:
	std::unordered_map<std::string, DNAread> myMap;

	hash_table();
	std::string get_key (std::string);
	bool has_key(std::string);
	bool add(std::string, std::string, std::string);
	std::string get_info (std::string, std::string, std::string, int);
	std::string get_seq(std::string);
	std::string get_qual(std::string);
	bool erase(std::string);
	bool is_empty();
	int size();
	void clear();
	~hash_table();
};

#endif