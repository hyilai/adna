# Final names of binaries
EXECUTABLE = Bin/adna
SO_LIBRARY = Bin/libzip.so

# C & C++ compiler
#CC       = gcc
#CXX      = g++-4.8
# CC        = clang
# CXX       = clang++
CC		  = /usr/local/gcc-5.3.0/bin/gcc
CXX 	  = /usr/local/gcc-5.3.0/bin/g++
CFLAGS    = -fPIC -Wno-enum-conversion -O3
CXXFLAGS  = -fPIC -std=c++11 -O3
MPICXX	  = mpic++ -std=c++11
SPMPICXX  = /usr/local/gcc-5.3.0/openmpi/bin/mpic++
MPIFLAG	  = -std=c++11

# Linker flags
LDFLAGS   = -pthread

# Sources of external libraries
SRC_ZLIB  = $(wildcard Source/ZipLib/extlibs/zlib/*.c)
SRC_LZMA  = $(wildcard Source/ZipLib/extlibs/lzma/unix/*.c)
SRC_BZIP2 = $(wildcard Source/ZipLib/extlibs/bzip2/*.c)

# ZipLib sources
SRC = \
		$(wildcard Source/ZipLib/*.cpp)        \
		$(wildcard Source/ZipLib/detail/*.cpp)

# Object files			
OBJS = \
		$(SRC:.cpp=.o)	   \
		$(SRC_ZLIB:.c=.o)  \
		$(SRC_LZMA:.c=.o)  \
		$(SRC_BZIP2:.c=.o)


# Rules
all: $(EXECUTABLE) $(SO_LIBRARY)

$(EXECUTABLE): $(OBJS)
	$(SPMPICXX) $(MPIFLAG) Source/SmithWaterman.cpp Source/global.cpp Source/hash_table.cpp Source/steps.cpp Source/memory_usage.cpp Source/utilities.cpp Source/MPI_readFastq2.cpp Source/main.cpp -o $@ $^

$(SO_LIBRARY): $(OBJS)
	$(CXX) $(LDFLAGS) -shared -o $@ $^

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

%.o: %.c
	$(CC) $(CFLAGS) -c $< -o $@

clean:
	rm -rf `find Source -name '*.o'` ziplib.tar.gz Bin/*.zip Bin/out* $(EXECUTABLE) $(SO_LIBRARY)

tarball:
	tar -zcvf ziplib.tar.gz *
	
