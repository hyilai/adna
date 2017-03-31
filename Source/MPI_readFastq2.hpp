#ifndef MPI_READ
#define MPI_READ

struct DNAread;


std::string MPI_receive_string (int, int);

void MPI_trim_file(char*, char*, int, int, int, std::vector<std::string>, bool);
void MPI_trim_and_match(char*, char*, int, int, int, std::vector<std::string>, bool);
double MPI_process_reads(char*, char*, char*, int, int, std::vector<std::string>, int, bool);
void trim_file(char*, char*, int, std::vector<std::string>, bool);
double process_reads (char*, char*, char*, int, std::vector<std::string>, int, bool);


#endif