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
	unordered_map<std::string, DNAread> temp;
	myMap = temp;
}


// get key for hash table
string hash_table::get_key (string info) {

	// tokenize
	stringstream ss;
	ss.str(info);
	string token;
	for (int i = 0; i < 5; i++) {
		getline(ss, token, ':');
	}
		
	string x;
	getline(ss, x, ':');
	string y;
	getline(ss, y, ' ');
	
	string key = x + ":" + y;

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

bool hash_table::add (string key, string info, string seq, string qual) {
	if (!has_key(key)) {
		DNAread temp;
		temp.info = info;
		temp.sequence = seq;
		temp.quality = qual;
		myMap.insert({key, temp});
		return true;
	} else {
		return false;
	}
}


string hash_table::get_info (string key) {
	return myMap[key].info;
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