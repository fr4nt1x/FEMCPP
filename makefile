CC=g++
CFLAGS=-c -Wall --std=c++11
LDFLAGS=-lCGAL -lgmp -larmadillo
SOURCES=run.cpp FiniteElementSolver.cpp
OBJECTS=$(SOURCES:.cpp=.o)
EXECUTABLE=./bin/run.exe

all: $(SOURCES) $(EXECUTABLE)
    
$(EXECUTABLE): $(OBJECTS) 
	$(CC) $(LDFLAGS) $(OBJECTS) -o $@

.cpp.o:
	$(CC) $(CFLAGS) $< -o $@

precompile:
	$(CC) $(CFLAGS) $(LDFLAGS) FiniteElementSolver.h -o FiniteElementSolver.h.gch 

clean:
	rm *o ./bin/run.exe	
