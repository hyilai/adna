CC = g++
SWH = SmithWaterman.hpp
SWCPP = SmithWaterman.cpp
TESTCPP = test.cpp
DNAGENCPP = DNAGeneratorV5.cpp

SW = SmithWaterman
TEST = test
DNAGEN = dnagen

all: $(SW) $(TEST) 

$(TEST): $(TESTCPP) $(SWCPP) $(SWH)
	$(CC) -g $(TESTCPP) $(SW).o -o $(TEST)

$(SW): $(SWCPP) $(SWH)
	$(CC) -c $(SWCPP) -o $(SW).o

$(DNAGEN): $(DNAGENCPP)
	$(CC) -g $(DNAGENCPP) -o $(DNAGEN)

