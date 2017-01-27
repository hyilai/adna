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
#define MISMATCH -15

//gap-scoring scheme
#define GAP_EXTENSION 10
#define GAP_PENALTY 20


class functions {
public:
	/**
	 **	this function calculates the individual scores
	 ** for each cell in the table
	 **/
	Grid build_grid (string str1, string str2) {

		int m = str1.length();
		int n = str2.length();
		grid.resize(m+1, vector<int>(n+1, 0));

		// cout << str1 << " " << str2 << " " << n << endl;

		int highest = 0;
		highest_i = 0;
		highest_j = 0;

		for (int i = 1; i < grid.size(); i++) {
			char a = str1[i-1];
			for (int j = 1; j < grid[i].size(); j++) {
				char b = str2[j-1];

				int score = 0;
				int match = 0;
				int deletion = 0;
				int insertion = 0;

				//match/mismatch
				int s = similarity(str1, str2, i, j);
				match = grid[i-1][j-1] + s;

				//deletion
				for (int k = i-1; k > 0; k--) {
					int temp = grid[i-k][j] + gap(k);
					if (deletion < temp) deletion = temp;
				}

				//insertion
				for (int l = j-1; l > 0; l--) {
					int temp = grid[i][j-l] + gap(l);
					if (deletion < temp) insertion = temp;
				}

				//find maximum between the three
				score = (score < match) ? match : score;
				score = (score < deletion) ? deletion : score;
				score = (score < insertion) ? insertion : score;

				grid[i][j] = score;

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
	
	int gap (int i) {
		return -(GAP_PENALTY + (GAP_EXTENSION * i) );
	}

	int similarity (string str1, string str2, int i, int j) {
		char a = str1[i-1];
		char b = str2[j-1];
		
		if (grid[i-1][j-1] > 0) {
			return (a == b) ? MATCH : MISMATCH;
		} else {
			return (a == b) ? MATCH - (j-1) : MISMATCH;
		}
	}
};


SmithWaterman::SmithWaterman (string str1, string str2, string q1, string q2, int match_score) : str1(str1), str2(str2), q1(q1), q2(q2), match_score(match_score) {
	
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

	trimmed = str1.substr(0, highest_i);
	trimmed_q = q1.substr(0, highest_i);
	matched = str1.substr(0, highest_i);
	string r = str1.substr(highest_i, str1.length()-highest_i);
	q1 = q1.substr(highest_i, q1.length()-highest_i);
	return r;
}


/**
 ** This function trims ending adapter from 
 ** sequence read
 **/	
string SmithWaterman::trim_from_ending () {

	//cout << highest_i << "," << highest_j << endl;
	//functions f;
	//f.print_grid(str1, str2, grid);

	//get trimmed sequence
	int curr_score = grid[highest_i][highest_j];	//highest score in the grid
	int i = highest_i;
	int j = highest_j;
	while (curr_score != 0 && i > 0 && j > 0) {
		int current = grid[i][j];
		int left = grid[i][j-1];
		int upper_left = grid[i-1][j-1];
		int upper = grid[i-1][j];

		int next_i = i-1;
		int next_j = j-1;
		int next_highest = 0;

		if (left > next_highest) {
			//left is biggest
			next_highest = left;
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
		}

		curr_score = next_highest;
		i = next_i;
		j = next_j;

	}

	trimmed = str1.substr(i,str1.length()-i);
	trimmed_q = q1.substr(i,q1.length()-i);
	matched = str1.substr(i,highest_i-i);
	
	string r = str1.substr(0,i);
	q1 = q1.substr(0,i);

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


/**
 ** WIP
 **/
bool SmithWaterman::match_reads () {
	int highest =  grid[highest_i][highest_j];

	// get matched sequence
	int curr_score = highest;
	int i = highest_i;
	int j = highest_j;
	while (curr_score != 0 && i > 0 && j > 0) {
		int current = grid[i][j];
		int left = grid[i][j-1];
		int upper_left = grid[i-1][j-1];
		int upper = grid[i-1][j];

		int next_i = i-1;
		int next_j = j-1;
		int next_highest = 0;

		if (left > next_highest) {
			//left is biggest
			next_highest = left;
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
		}

		curr_score = next_highest;
		i = next_i;
		j = next_j;

	}

	if (highest_j - j < match_score || highest_i - i < match_score) {
		return false;
	} else {
		// concatenate matching strings
		stringstream ss;
		if (j < i) {
			ss << str1.substr(0,highest_i);
			ss << str2.substr(highest_j, str2.length() - highest_j);
		} else {
			ss << str2.substr(0,highest_j);
			ss << str1.substr(highest_i, str1.length() - highest_i);
		}
		matched = ss.str();

		return true;
	}
}

// SmithWaterman::~SmithWaterman() {
// }

// extern "C" SmithWaterman* create(string str1, string str2, int match_score) {
// 	return new SmithWaterman(str1, str2, match_score);
// }


// extern "C" void erase(SmithWaterman* sw) {
// 	delete sw;
// }
