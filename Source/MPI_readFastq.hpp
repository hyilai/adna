#ifndef MPI_READ
#define MPI_READ


void MPI_trim_file(char*, char*, int, int, std::vector<std::string>, bool);
void MPI_trim_and_match(char*, char*, int, int, std::vector<std::string>, bool);
void MPI_process_reads(char*, char*, char*, int, int, std::vector<std::string>, bool);
void trim_file(char*, char*, std::vector<std::string>, bool);
void process_reads (char*, char*, char*, std::vector<std::string>, bool);


#endif