# Final names of binaries
EXECUTABLE = Bin/adna

# C & C++ compiler
#CC       = gcc
#CXX      = g++-4.8
# CC        = clang
# CXX       = clang++
CC		  = /usr/local/gcc-5.3.0/bin/gcc
CXX 	  = /usr/local/gcc-5.3.0/bin/g++
CXXFLAGS  = -std=c++11
# MPICXX	  = mpic++
MPICXX    = /usr/local/gcc-5.3.0/openmpi/bin/mpic++
MPIFLAG	  = -std=c++11
BOOSTFLAG = -lz -lboost_iostreams

$SWCPP 	  = Source/SmithWaterman.cpp
$SWH	  = Source/SmithWaterman.hpp

$HASHCPP  = Source/hash_table.cpp
$HASHH 	  = Source/hash_table.hpp


# Rules
all: $(EXECUTABLE)

$(EXECUTABLE):
	$(MPICXX) $(MPIFLAG) Source/SmithWaterman.cpp Source/global.cpp Source/hash_table.cpp Source/steps.cpp Source/memory_usage.cpp Source/utilities.cpp Source/MPI_readFastq.cpp Source/main.cpp -o $@ $^


clean:


tarball:
	tar -zcvf ziplib.tar.gz *
	
