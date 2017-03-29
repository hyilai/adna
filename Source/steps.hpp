#ifndef H_SW
#define H_SW

#include "SmithWaterman.hpp"
class SmithWaterman;


bool strip_t(std::string&, std::string&);
std::string trim_read(std::string, std::string, std::string&, std::string&, std::string&, std::vector<std::string>);
// std::string trim_read(std::string, std::string, std::string&, std::string&, std::string&, std::vector<std::string>);
std::string concat_reads(std::string, std::string, std::string, std::string, std::string&);
bool quality_check(std::string);

#endif