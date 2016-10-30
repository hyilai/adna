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
	SmithWaterman (std::string, std::string, int = 0);
	std::string trim_both_sides ();
	std::string trim_from_beginning ();
	std::string trim_from_ending ();
	// ~SmithWaterman();
private:
	std::string str1;
	std::string str2;
	std::vector<std::vector<int> > grid;
	int match_score;
};

#endif