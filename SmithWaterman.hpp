#ifndef SW_HPP
#define SW_HPP

#include <cstdlib>
#include <cstdio>
#include <string.h>
#include <vector>
#include <iomanip>
#include <iostream>
// #include <pthread.h>


class SmithWaterman {
public:
	SmithWaterman (std::string, std::string, std::string, std::string, int = 0);
	int get_highest_score();
	std::string trim_both_sides ();
	std::string trim_from_beginning ();
	std::string trim_from_ending ();
	std::string get_quality1();
	std::string get_quality2();
	std::string get_trimmed();
	std::string get_matched_string();
	std::string match_reads();
	// ~SmithWaterman();
private:
	std::string str1;
	std::string str2;
	std::string trimmed;
	std::string matched;
	std::string q1;		//quality string for str1
	std::string q2;		//quality string for str2, would be empty if comparing sequence to adapter
	int match_score;	//minimum score for matching two reads
	std::vector<std::vector<int> > grid;
};

#endif