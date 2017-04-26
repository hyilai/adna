/************************************************************
 *
 *	This portion uses the Smith Waterman Algorithm 
 *	to do several tasks:
 *		1. find adapter fragments within the reads
 *			and trim them out
 *		2. find matches between two sequences by
 *			using their quality control strings
 *
 ************************************************************/


#include "SmithWaterman.hpp"

using namespace std;

typedef int* Grid;

//similarity function
#define MATCH 5
#define MISMATCH -8

//gap-scoring scheme
#define GAP_EXTENSION 10
#define GAP_PENALTY 10


SmithWaterman::SmithWaterman(int match_length): match_length(match_length) {
	str1 = "";
	str2 = "";
	q1 = "";
	q2 = "";
	trimmed = "";
	trimmed_q = "";
	concat = "";
	concat_q = "";
	highest_i = 0;
	highest_j = 0;
	m = 0;
	n = 0;
}


SmithWaterman::SmithWaterman(int max_rows, int max_cols, int match_length): match_length(match_length) {
	set_grid_size(max_rows, max_cols);

	str1 = "";
	str2 = "";
	q1 = "";
	q2 = "";
	trimmed = "";
	trimmed_q = "";
	concat = "";
	concat_q = "";
	highest_i = 0;
	highest_j = 0;
	m = 0;
	n = 0;
}


bool SmithWaterman::set_grid_size(int rows, int cols) {
	if (rows == 0 || cols == 0) {
		cout << "Error: grid size should be at least 1 x 1" << endl;
		return false;
	}

	// make a 2D grid masked as a 1D array
	grid = new int[(rows+1) * (cols+1)];
	return true;
}

bool SmithWaterman::set_strings (string s1, string qual1, string s2) {
	str1 = s1;
	str2 = s2;
	m = s1.length();
	n = s2.length();
	q1 = qual1;
	q2 = "";
	return m > 0 && n > 0;
}

bool SmithWaterman::set_strings (string s1, string qual1, string s2, string qual2) {
	str1 = s1;
	str2 = s2;
	m = s1.length();
	n = s2.length();
	q1 = qual1;
	q2 = qual2;
	return m > 0 && n > 0;
}


/**
 **	this function calculates the individual scores
 ** for each cell in the table
 **/
void SmithWaterman::build_grid () {

	int highest = 0;
	highest_i = 0;
	highest_j = 0;

	// for each cell in the grid
	for (int i = 0; i < m+1; i++) {

		char a;
		if (i > 0) a = str1[i-1];

		for (int j = 0; j < n+1; j++) {

			if (i == 0 || j == 0) {
				grid[get_index(i, j)] = 0;
				continue;
			}

			char b = str2[j-1];

			int score = 0;

			int top = grid[get_index(i-1, j)];
			int top_left = grid[get_index(i-1, j-1)];
			int left = grid[get_index(i, j-1)];

			if (top == 0 && top_left == 0 && left == 0) {

				grid[get_index(i, j)] = max(0, similarity(i, j, top_left));

			} else {

				int match = 0, deletion = 0, insertion = 0;

				// match/mismatch score
				match = top_left + similarity(i, j, top_left);

				// deletion score
				deletion = left - GAP_PENALTY;

				// insertion score
				insertion = top - GAP_PENALTY;

				// find maximum between the three; minimum score is 0
				score = (score < match) ? match : score;
				score = (score < deletion) ? deletion : score;
				score = (score < insertion) ? insertion : score;

				// set cell to the maximum score
				grid[get_index(i, j)] = score;
			}

			// remember the highest score in the grid
			if (score > highest) {
				highest = score;
				highest_i = i;
				highest_j = j;
			}
		}
	}
}

void SmithWaterman::print_grid () {
	cout << "x" << "\t" << setw(6) << "- ";
	for (int i = 0; i < str2.length(); i++) {
		cout << setw(3) << str2[i] << " ";
	}
	cout << endl;
	for (int i = 0; i <= m; i++) {
		if (i != 0) {
			cout << i << "\t" << str1[i-1] << " ";
		} else {
			cout << i << "\t" << "- ";
		}
		for (int j = 0; j <= n; j++) {
			cout << setw(3) << grid[get_index(i, j)] << " ";
		}
		cout << endl;
	}
	cout << endl;
}


int SmithWaterman::get_highest_i () {
	return highest_i;
}

int SmithWaterman::get_highest_j () {
	return highest_j;
}


int SmithWaterman::gap (int i, int prev_score) {
	// int n = prev_score / MATCH;
	// return (GAP_PENALTY / (n + 1) + GAP_EXTENSION) * i;
	return GAP_PENALTY + GAP_EXTENSION * i;
}

int SmithWaterman::similarity (int i, int j, int prev) {
	char a = str1[i-1];
	char b = str2[j-1];
	
	// if we're matching against the first character in str2
	if (j - 1 == 0) {
		return (a == b) ? MATCH : MISMATCH;
	} else {

		// if the last set of characters matched between str1 and str2
		if (prev > 0) {
			return (a == b) ? MATCH : MISMATCH;
		} else {	// else return a harsher score
			return (a == b) ? MATCH - j*j : MISMATCH;
		}
	}
}

int SmithWaterman::get_index (int r, int c) {
	return r*(n+1) + c;
}

int SmithWaterman::get_highest_score () {
	return grid[get_index(highest_i, highest_j)];
}

/**
 ** This function trims beginning adapter from 
 ** sequence read
 **/	
string SmithWaterman::trim_from_beginning () {

	// string r = str1;
	if (highest_i >= match_length) {
		trimmed = str1.substr(0, highest_i);
		trimmed_q = q1.substr(0, highest_i);
		str1 = str1.substr(highest_i, str1.length()-highest_i);
		q1 = q1.substr(highest_i, q1.length()-highest_i);
	} else {
		trimmed = "";
		trimmed_q = "";
	}
	return str1;
}


/**
 ** This function trims ending adapter from sequence read
 **/	
string SmithWaterman::trim_from_ending () {

	int curr_score = grid[get_index(highest_i, highest_j)];	//highest score in the grid
	int i = highest_i;
	int j = highest_j;
	while (curr_score != 0 && i > 0 && j > 0) {
		int current = grid[get_index(i, j)];
		int left = grid[get_index(i, j-1)];
		int upper_left = grid[get_index(i-1, j-1)];
		int upper = grid[get_index(i-1, j)];

		int next_i = i;
		int next_j = j;
		int next_highest = 0;

		if (left > next_highest) {
			//left is biggest
			next_highest = left;
			next_i = i;
			next_j = j-1;
		}
		if (upper_left > next_highest) {
			//upper_left is biggest
			next_highest = upper_left;
			next_i = i-1;
			next_j = j-1;
		}
		if (upper > next_highest) {
			//upper is bigest
			next_highest = upper;
			next_i = i-1;
			next_j = j;
		}

		curr_score = next_highest;
		i = next_i;
		j = next_j;
	}

	// string r = str1;

	// trim only if the length of the matching substring is more than the threshold
	// everything after the matching substring is thrown out as well
	if (highest_i - i >= match_length) {
		trimmed = str1.substr(i-1, str1.length() - i);
		trimmed_q = q1.substr(i-1, q1.length() - i);
		
		str1 = str1.substr(0, i);
		q1 = q1.substr(0, i);
	} else {
		trimmed = "";
		trimmed_q = "";
	}

	return str1;
}

string SmithWaterman::get_quality1 () {
	return q1;
}

string SmithWaterman::get_trimmed_quality () {
	return trimmed_q;
}

string SmithWaterman::get_quality2 () {
	return q2;
}

string SmithWaterman::get_trimmed() {
	return trimmed;
}

string SmithWaterman::get_concat_string() {
	return concat;
}

string SmithWaterman::get_concat_quality() {
	return concat_q;
}


/**
 ** find the matching substring between the two reads and concatenate them
 **/
bool SmithWaterman::concat_strings () {

	int curr_score = grid[get_index(highest_i, highest_j)];	// highest score in the grid
	int i = highest_i;
	int j = highest_j;
	int prev_i = i;
	int prev_j = j;

	// go back through the grid to find the matching subsequence
	while (curr_score != 0 && i > 0 && j > 0) {
		int current = grid[get_index(i, j)];
		int left = grid[get_index(i, j-1)];
		int upper_left = grid[get_index(i-1, j-1)];
		int upper = grid[get_index(i-1, j)];

		int next_i = i;
		int next_j = j;
		int next_highest = 0;

		if (left > next_highest) {
			//left is biggest
			next_highest = left;
			next_i = i;
			next_j = j-1;
		}
		if (upper_left > next_highest) {
			//upper_left is biggest
			next_highest = upper_left;
			next_i = i-1;
			next_j = j-1;
		}
		if (upper > next_highest) {
			//upper is biggest
			next_highest = upper;
			next_i = i-1;
			next_j = j;
		}

		curr_score = next_highest;
		prev_i = i;
		prev_j = j;
		i = next_i;
		j = next_j;
	}


	// if the length of the matching subsequence is less than the minimum match length do nothing
	if (highest_j - (prev_j - 1) < match_length || highest_i - (prev_i - 1) < match_length) {
		return false;
	} else {	// else concatenate the two strings togther
		// concatenate matching strings
		stringstream ss;
		stringstream ss_q;
		ss << str1.substr(0, i-1);
		ss << str2.substr(0, str2.length());
		ss_q << q1.substr(0, i-1);
		ss_q << q2.substr(0, str2.length());
		
		concat = ss.str();
		concat_q = ss_q.str();

		return true;
	}
}

// destructor
SmithWaterman::~SmithWaterman() {
	delete [] grid;
}
