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

typedef vector<vector<int> > Grid;

//similarity function
#define MATCH 5
#define MISMATCH -8

//gap-scoring scheme
#define GAP_EXTENSION 10
#define GAP_PENALTY 7


class functions {
public:
	/**
	 **	this function calculates the individual scores
	 ** for each cell in the table
	 **/
	Grid build_grid (string str1, string str2) {

		int m = str1.length();
		int n = str2.length();

		// make a grid of size m by n
		grid.resize(m+1, vector<int>(n+1, 0));

		int highest = 0;
		highest_i = 0;
		highest_j = 0;

		// for each cell in the grid
		for (int i = 1; i < grid.size(); i++) {
			char a = str1[i-1];
			for (int j = 1; j < grid[i].size(); j++) {
				char b = str2[j-1];

				int score = 0, match = 0, deletion = 0, insertion = 0;

				// match/mismatch score
				int s = similarity(str1, str2, i, j);
				match = grid[i-1][j-1] + s;

				// deletion score
				for (int k = i-1; k > 0; k--) {
					int temp = grid[i-k][j] - gap(k, grid[i-k][j]);
					if (deletion < temp) deletion = temp;
				}

				// insertion score
				for (int l = j-1; l > 0; l--) {
					int temp = grid[i][j-l] - gap(l, grid[i][j-l]);
					if (insertion < temp) insertion = temp;
				}

				// find maximum between the three; minimum score is 0
				score = (score < match) ? match : score;
				score = (score < deletion) ? deletion : score;
				score = (score < insertion) ? insertion : score;

				// set cell to the maximum score
				grid[i][j] = score;

				// remember the highest score in the grid
				if (score > highest) {
					highest = score;
					highest_i = i;
					highest_j = j;
				}
			}
		}
		return grid;
	}

	void print_grid (string str1, string str2, Grid grid) {
		int m = str1.length();
		int n = str2.length();
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
				cout << setw(3) << grid[i][j] << " ";
			}
			cout << endl;
		}
		cout << endl;
	}

	int get_highest_i () {
		return highest_i;
	}

	int get_highest_j () {
		return highest_j;
	}

private:
	Grid grid;
	int highest_i;
	int highest_j;
	
	int gap (int i, int prev_score) {
		int n = prev_score / MATCH;
		// original: GAP_PENALTY + GAP_EXTENSION * i
		return (GAP_PENALTY / (n + 1) + GAP_EXTENSION) * i;
	}

	int similarity (string str1, string str2, int i, int j) {
		char a = str1[i-1];
		char b = str2[j-1];
		
		// if we're matching against the first character in str2
		if (j - 1 == 0) {
			return (a == b) ? MATCH : MISMATCH;
		} else {

			// if the last set of characters matched between str1 and str2
			if (grid[i-1][j-1] > 0) {
				return (a == b) ? MATCH : MISMATCH;
			} else {	// else return a harsher score
				return (a == b) ? MATCH - j*j : MISMATCH;
			}
		}
	}
};


SmithWaterman::SmithWaterman (string str1, string str2, string q1, int match_length) : str1(str1), str2(str2), q1(q1), match_length(match_length) {
	
	q2 = "";
	trimmed = "";
	functions f;

	int m = str1.length();
	int n = str2.length();

	if (m == 0 || n == 0) {
		cout << "Error: string(s) length is zero" << endl;
		exit(EXIT_FAILURE);
	}

	grid = f.build_grid(str1, str2);

	highest_i = f.get_highest_i();
	highest_j = f.get_highest_j();

	//f.print_grid(str1,str2,grid);
}


SmithWaterman::SmithWaterman (string str1, string str2, string q1, string q2, int match_length) : str1(str1), str2(str2), q1(q1), q2(q2), match_length(match_length) {
	
	trimmed = "";
	functions f;

	int m = str1.length();
	int n = str2.length();

	if (m == 0 || n == 0) {
		cout << "Error: string(s) length is zero" << endl;
		exit(EXIT_FAILURE);
	}

	grid = f.build_grid(str1, str2);

	highest_i = f.get_highest_i();
	highest_j = f.get_highest_j();

	//f.print_grid(str1,str2,grid);
}

int SmithWaterman::get_highest_score () {
	return grid[highest_i][highest_j];
}

string SmithWaterman::trim_both_sides () {
	string trimmed = trim_from_beginning();
	str1 = trimmed;
	functions f;
	grid = f.build_grid(str1, str2);

	//f.print_grid(str1,str2,grid);

	return trim_from_ending();
}

/**
 ** This function trims beginning adapter from 
 ** sequence read
 **/	
string SmithWaterman::trim_from_beginning () {

	string r = str1;
	if (highest_i >= match_length) {
		trimmed = str1.substr(0, highest_i);
		trimmed_q = q1.substr(0, highest_i);
		matched = str1.substr(0, highest_i);
		r = str1.substr(highest_i, str1.length()-highest_i);
		q1 = q1.substr(highest_i, q1.length()-highest_i);
	} else {
		trimmed = "";
		trimmed_q = "";
		matched = "";
	}
	return r;
}


/**
 ** This function trims ending adapter from sequence read
 **/	
string SmithWaterman::trim_from_ending () {

	//get trimmed sequence
	int curr_score = grid[highest_i][highest_j];	//highest score in the grid
	int i = highest_i;
	int j = highest_j;
	while (curr_score != 0 && i > 0 && j > 0) {
		int current = grid[i][j];
		int left = grid[i][j-1];
		int upper_left = grid[i-1][j-1];
		int upper = grid[i-1][j];

		int next_i, next_j;
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

	string r = str1;

	// trim only if the length of the matching substring is more than the threshold
	// everything after the matching substring is thrown out as well
	if (highest_i - i >= match_length) {
		trimmed = str1.substr(i-1, str1.length() - i);
		trimmed_q = q1.substr(i-1, q1.length() - i);
		matched = str1.substr(i-1, highest_i - i);
		
		r = str1.substr(0, i);
		q1 = q1.substr(0, i);
	} else {
		trimmed = "";
		trimmed_q = "";
		matched = "";
	}

	// cout << i << " " << r << endl;
	return r;
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

string SmithWaterman::get_matched_string() {
	return matched;
}

string SmithWaterman::get_matched_quality() {
	return matched_q;
}


/**
 ** find the matching substring between the two reads and concatenate them
 **/
bool SmithWaterman::match_reads () {

	// get the highest score in the grid
	int highest =  grid[highest_i][highest_j];

	// go back through the grid to find the matching subsequence
	int curr_score = highest;
	int i = highest_i;
	int j = highest_j;
	int prev_i = i;
	int prev_j = j;
	while (curr_score != 0 && i > 0 && j > 0) {
		int current = grid[i][j];
		int left = grid[i][j-1];
		int upper_left = grid[i-1][j-1];
		int upper = grid[i-1][j];

		int next_i, next_j;
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
		
		matched = ss.str();
		matched_q = ss_q.str();

		return true;
	}
}

// destructor
SmithWaterman::~SmithWaterman() {

}
