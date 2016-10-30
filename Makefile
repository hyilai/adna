CC = g++
FILE1 = DNAGeneratorV3.cpp
FILE2 = randomStringGeneratorV1.cpp

PROG1 = dnaGen
PROG2 = randStr


all: $(PROG1) $(PROG2)

$(PROG1): $(FILE1)
	$(CC) $(FILE1) -o $(PROG1)

$(PROG2): $(FILE2)
	$(CC) $(FILE2) -o $(PROG2)

