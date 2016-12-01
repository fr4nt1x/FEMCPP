CC=g++
CFLAGS=-c -Wall --std=c++11
LDFLAGS=-lCGAL -lgmp -larmadillo

SOURCES=run.cpp FiniteElementSolver.cpp
SOURCES1=run1.cpp FiniteElementSolver.cpp

OBJECTS=$(SOURCES:.cpp=.o)
OBJECTS1=$(SOURCES1:.cpp=.o)
EXECUTABLE=./bin/run.exe
EXECUTABLE1=./bin/run1.exe

all: $(SOURCES) $(EXECUTABLE)
    
all1: $(SOURCES1) $(EXECUTABLE1)

$(EXECUTABLE): $(OBJECTS) 
	$(CC) $(LDFLAGS) $(OBJECTS) -o $@

$(EXECUTABLE1): $(OBJECTS1) 
	$(CC) $(LDFLAGS) $(OBJECTS1) -o $@

.cpp.o:
	$(CC) $(CFLAGS) $< -o $@

precompile:
	$(CC) $(CFLAGS) $(LDFLAGS) FiniteElementSolver.h -o FiniteElementSolver.h.gch 

clean:
	rm *o ./bin/run.exe	
