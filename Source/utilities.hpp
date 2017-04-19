#ifndef H_UTIL
#define H_UTIL

void merge_files (std::string,  int);
void merge_worker_files (std::string, int);

std::string get_file_name(std::string, int, std::string);
std::string get_file_name(std::string, int, int, std::string);

bool is_same_pair(std::string, std::string);

void write_to_fastq(std::ofstream&, std::string, std::string, std::string, std::string);

void print_diagnostics(double, int, int, int, int, int, int);

#endif