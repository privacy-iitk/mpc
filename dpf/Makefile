CXX = g++
CXXFLAGS = -std=c++17 -march=native -O3 -g -Wall -pedantic -Wno-ignored-attributes -fopenmp -msse2 -mssse3 -DBOOST_ERROR_CODE_HEADER_ONLY
LIBS = -lbsd -fopenmp -lboost_system -lboost_thread -pthread

all: dpf

dpf: dpf.cpp dpf.h prg.h block.h
	$(CXX) $(CXXFLAGS) -o $@ dpf.cpp $(LIBS)

preproc: preproc.cpp dpf.h prg.h block.h
	$(CXX) $(CXXFLAGS) -o $@ preproc.cpp $(LIBS)
clean:
	rm -f dpf
