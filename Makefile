# Final names of binaries
EXECUTABLE = Bin/adna
RUNTIME	   = Bin/test_runtime

# C & C++ compiler
CC		  = /usr/local/gcc-5.3.0/bin/gcc
CXX 	  = /usr/local/gcc-5.3.0/bin/g++
CXXFLAGS  = -std=c++11
# MPICXX	  = mpic++
MPICXX    = /usr/local/gcc-5.3.0/openmpi/bin/mpic++
MPIFLAG	  = -std=c++11
BOOSTFLAG = -lz -lboost_iostreams

$SWCPP 	  = Source/SmithWaterman.cpp
$SWH	  = Source/SmithWaterman.hpp


# Rules
all: $(EXECUTABLE)

$(EXECUTABLE):
	$(MPICXX) $(MPIFLAG) Source/SmithWaterman.cpp Source/global.cpp Source/steps.cpp Source/memory_usage.cpp Source/utilities.cpp Source/MPI_readFastq.cpp Source/main.cpp -o $@ $^

$(RUNTIME):
	$(MPICXX) $(MPIFLAG) Source/SmithWaterman.cpp Source/global.cpp Source/steps.cpp Source/memory_usage.cpp Source/utilities.cpp Source/MPI_readFastq2.cpp Source/main2.cpp -o $@ $^

clean:
	rm $(EXECUTABLE) $(RUNTIME)
