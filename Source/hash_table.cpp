#include <string.h>
#include <sstream>
#include <fstream>
#include <iostream>
#include <unordered_map>
#include <vector>
#include "hash_table.hpp"

using namespace std;


// constructor
hash_table::hash_table() {

}


// get key for hash table
string hash_table::get_key (string info) {

	// tokenize
	string flowcell, x, y;

	size_t found = info.find(" ");
	info = info.substr(0, found);

	found = info.find_last_of(":");			// y coordinate
	y = info.substr(found+1);
	info = info.substr(0, found);
	found = info.find_last_of(":");			// x coordinate
	x = info.substr(found+1);
	info = info.substr(0, found);
	found = info.find_last_of(":");			// flowcell
	flowcell = info.substr(found+1);
	info = info.substr(0, found);
	
	string key = flowcell + ":" + x + ":" + y;

	return key;	
}

bool hash_table::has_key (string key) {
	unordered_map<string, DNAread>::const_iterator found = myMap.find(key);
	if (found == myMap.end()) {
		return false;
	} else {
		return true;
	}
}

bool hash_table::add (string key, string seq, string qual) {
	if (!has_key(key)) {
		DNAread temp;
		temp.sequence = seq;
		temp.quality = qual;
		myMap.insert({key, temp});
		return true;
	} else {
		return false;
	}
}

string hash_table::get_info (string first, string second, string key, int num) {
	stringstream ss;
	ss << first << ":" << key << " " << num << ":" << second;
	return ss.str();
}


string hash_table::get_seq (string key) {
	return myMap[key].sequence;
}

string hash_table::get_qual (string key) {
	return myMap[key].quality;
}

bool hash_table::erase(string key) {
	myMap.erase(key);
}

bool hash_table::is_empty() {
	return myMap.empty();
}

int hash_table::size() {
	return myMap.size();
}

void hash_table::clear(){
	myMap.clear();
}


// destructor
hash_table::~hash_table() {
	clear();
}