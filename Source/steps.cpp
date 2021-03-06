#include <cstdlib>
#include <cstdio>
#include <string.h>
#include <string>
#include <vector>
#include "global.hpp"
#include "steps.hpp"

using namespace std;


// step 1
// change quality score to 2 ('#') if the leading and/or last 2 bases are 'T's
bool strip_t (string &sequence, string &quality) {
	bool has_t = false;
	int seq_length = sequence.length();
	if (sequence[0] == 'T') {
		quality[0] = '#';
		has_t = true;
	}
	if (sequence[seq_length-2] == 'T' && sequence[seq_length-1] == 'T') {
		quality[seq_length-2] = quality[seq_length-1] = '#';
		has_t = true;
	}
	return has_t;
}

// step 2 if no adapter file is specified
string trim_read (SmithWaterman *sw, string sequence, string quality, string &new_quality, string &junk, string &junk_quality, vector<string> adapters) {
	string read = sequence;
	int highest_score = 0;


	for (int i = 0; i < adapters.size(); i++) {

		sw->set_strings(sequence, quality, adapters[i]);
		sw->build_grid();

		int curr_high_score = sw->get_highest_score();

		if (highest_score < curr_high_score) {
			highest_score = curr_high_score;
			read = sw->trim_from_ending();
			new_quality = sw->get_quality1();
			junk = sw->get_trimmed();
			junk_quality = sw->get_trimmed_quality();
		}
	}


	return read;
}


// step 3
// merge two reads together
string concat_reads (SmithWaterman *sw, string sequence1, string sequence2, string quality1, string quality2, string &new_quality) {
	string concat = "";
	new_quality = "";

	sw->set_strings(sequence1, quality1, sequence2, quality2);
	sw->build_grid();
	// int curr_high_score = sw->get_highest_score();
	if (sw->concat_strings()) {
		concat = sw->get_concat_string();
		new_quality = sw->get_concat_quality();
	}
	
	return concat;
}


// step 4
// finds bases with score less than 15 (more than 5: fails)
bool quality_check (string quality) {
	int num_failed = 0;
	for (int i = 0; i < quality.length(); i++) {
		char cur = quality[i];
		if (int(cur) - 33 < 15) {	// adjust ascii value for base 33 quality strings
			num_failed++;
		}
		if (num_failed > 5) return false;
	}
	return true;
}
