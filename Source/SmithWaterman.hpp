#ifndef SW_HPP
#define SW_HPP

#include <cstdlib>
#include <cstdio>
#include <string.h>
#include <vector>
#include <iomanip>
#include <sstream>
#include <iostream>


class SmithWaterman {
public:
	// SmithWaterman (std::string, std::string, std::string, int = 0);
	// SmithWaterman (std::string, std::string, std::string, std::string, int = 0);
	SmithWaterman(int);
	SmithWaterman(int, int, int);

	bool set_grid_size(int, int);

	bool set_strings(std::string, std::string, std::string);
	bool set_strings(std::string, std::string, std::string, std::string);

	void build_grid();
	void print_grid();

	int get_highest_score();
	std::string get_quality1();
	std::string get_quality2();
	std::string get_trimmed();
	std::string get_trimmed_quality();
	std::string get_concat_string();
	std::string get_concat_quality();

	std::string trim_from_beginning ();
	std::string trim_from_ending ();
	bool concat_strings();
	~SmithWaterman();
private:
	std::string str1;	//sequence1
	std::string str2;	//sequence2
	std::string q1;		//quality string for str1
	std::string q2;		//quality string for str2, would be empty if comparing sequence to adapter
	std::string trimmed;	//sequence that is trimmed off
	std::string trimmed_q;	//quality string that is trimmed off
	std::string concat;	// the concatenated string
	std::string concat_q;
	int match_length;	//minimum matching length for two reads
	int* grid;			// score grid
	int highest_i;		// row number of the cell containing the highest score
	int highest_j;		// column of the cell containing the highest score
	int m;		// length of string 1
	int n;		// length of string 2

	int get_index(int, int);
	int get_highest_i();
	int get_highest_j();

	// functions to help building the grid
	int gap(int, int);
	int similarity(int, int, int);
};

#endif
