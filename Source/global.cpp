#include <string.h>
#include <string>
#include <fstream>
#include <vector>
#include "global.hpp"

using namespace std;

int minimum_read_length;
int minimum_match_length;

void set_read_length(int length) {
	minimum_read_length = length;
}

void set_match_length(int length) {
	minimum_match_length = length;
}

