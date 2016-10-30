#include "SmithWaterman.hpp"

using namespace std;

typedef vector<vector<int> >Grid;

//similarity function
#define MATCH 7
#define MISMATCH -4
// const int MATCH = 2;
// const int MISMATCH = -3;

//gap-scoring scheme
#define GAP_EXTENSION 8
#define GAP_PENALTY 10
// const int GAP = -4;


class functions {
public:
	/**
	 **	this function calculates the individual scores
	 ** for each cell in the table
	 **/
	Grid build_grid (string str1, string str2) {

		int m = str1.length();
		int n = str2.length();
		Grid grid(m+1, vector<int>(n+1, 0));

		// cout << str1 << " " << str2 << " " << n << endl;

		for (int i = 1; i < grid.size(); i++) {
			char a = str1[i-1];
			for (int j = 1; j < grid[i].size(); j++) {
				char b = str2[j-1];

				int score = 0;
				int match = 0;
				int deletion = 0;
				int insertion = 0;

				//match/mismatch
				int s = (a == b) ? MATCH : MISMATCH;
				match = grid[i-1][j-1] + s;

				//deletion
				for (int k = i-1; k > 0; k--) {
					int temp = grid[i-k][j] - gap(k);
					if (deletion < temp) deletion = temp;
				}
				// deletion = grid[i-1][j] + GAP;

				//insertion
				for (int l = j-1; l > 0; l--) {
					int temp = grid[i][j-l] - gap(l);
					if (deletion < temp) insertion = temp;
				}
				// insertion = grid[i][j-1] + GAP;

				//find maximum between the three
				score = (score < match) ? match : score;
				score = (score < deletion) ? deletion : score;
				score = (score < insertion) ? insertion : score;

				grid[i][j] = score;
			}
		}
		return grid;
	}

	void print_grid (string str1, string str2, Grid grid) {
		int m = str1.length();
		int n = str2.length();
		cout << setw(5) << "- ";
		for (int i = 0; i < str2.length(); i++) {
			cout << setw(2) << str2[i] << " ";
		}
		cout << endl;
		for (int i = 0; i < m+1; i++) {
			if (i != 0) {
				cout << str1[i-1] << " ";
			} else {
				cout << "- ";
			}
			for (int j = 0; j < n+1; j++) {
				cout << setw(2) << grid[i][j] << " ";
			}
			cout << endl;
		}
		cout << endl;
	}
private:
	int gap (int i) {
		return GAP_PENALTY + (GAP_EXTENSION * i);
	}

	// int similarity_function() {

	// }
};


SmithWaterman::SmithWaterman (string str1, string str2, int match_score) : str1(str1), str2(str2), match_score(match_score) {
	functions f;

	int m = str1.length();
	int n = str2.length();
	// grid.resize(m+1,vector<int>(n+1,0));
	if (m == 0 || n == 0) {
		cout << "Error: string(s) length is zero" << endl;
		exit(EXIT_FAILURE);
	}

	grid = f.build_grid(str1, str2);

	f.print_grid(str1,str2,grid);
}

string SmithWaterman::trim_both_sides () {
	string trimmed = trim_from_beginning();
	str1 = trimmed;
	functions f;
	grid = f.build_grid(str1, str2);

	int m = str1.length();
	int n = str2.length();
	cout << setw(5) << "- ";
	for (int i = 0; i < str2.length(); i++) {
		cout << setw(2) << str2[i] << " ";
	}
	cout << endl;
	for (int i = 0; i < m+1; i++) {
		if (i != 0) {
			cout << str1[i-1] << " ";
		} else {
			cout << "- ";
		}
		for (int j = 0; j < n+1; j++) {
			cout << setw(2) << grid[i][j] << " ";
		}
		cout << endl;
	}
	cout << endl;

	return trim_from_ending();
}

/**
 ** This function trims beginning adapter from 
 ** sequence read
 **/	
string SmithWaterman::trim_from_beginning () {
	int highest = 0;
	int highest_i = 0;
	int highest_j = 0;

	for (int i = 0; i < grid.size(); i++) {
		for (int j = 0; j <grid[i].size(); j++) {
			if (grid[i][j] > highest) {
				highest = grid[i][j];
				highest_i = i;
				highest_j = j;
			}
		}
	}

	string trimmed = str1.substr(highest_i, str1.length()-highest_i);
	return trimmed;
}

/**
 ** This function trims ending adapter from 
 ** sequence read
 **/	
string SmithWaterman::trim_from_ending () {
	int highest = 0;
	int highest_i = 0;
	int highest_j = 0;

	for (int i = 0; i < grid.size(); i++) {
		for (int j = 0; j <grid[i].size(); j++) {
			if (grid[i][j] >= highest) {
				highest = grid[i][j];
				highest_i = i;
				highest_j = j;
			}
		}
	}

	// cout << highest_i << "," << highest_j << endl;

	//get trimmed sequence
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
			// if (left != current) {
				//left is biggest
				next_highest = left;
				next_j = j-1;
			// }
		}
		if (upper_left > next_highest) {
			// if (upper_left != current) {
				//upper_left is biggest
				next_highest = upper_left;
				next_i = i-1;
				next_j = j-1;
			// }
		}
		if (upper > next_highest) {
			// if (upper != current) {
				//upper is bigest
				next_highest = upper;
				next_i = i-1;
			// }
		}

		i = next_i;
		j = next_j;

		// cout << i << "," << j << endl;
	}
	string r = str1.substr(0,i);

	// cout << i << "," << j << endl;

	// cout << i << " " << r << endl;
	return r;
}

/**
 ** 
 ** 
 **/
// string SmithWaterman::concatString () {

// }

// SmithWaterman::~SmithWaterman() {
// }

// extern "C" SmithWaterman* create(string str1, string str2, int match_score) {
// 	return new SmithWaterman(str1, str2, match_score);
// }


// extern "C" void erase(SmithWaterman* sw) {
// 	delete sw;
// }