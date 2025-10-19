#-*-Makefile-*- for ruferr
##
CXX=/opt/homebrew/bin/g++
CXXFLAGS=-std=c++20 -g -Wall -Wno-deprecated-declarations -O3 -fopenmp
CPPFLAGS=-I/opt/homebrew/include
LDFLAGS=-L/opt/homebrew/lib
LIBS=-lfftw3 -lfftw3_omp
all: ruferr

ruferr: ruferr.o coordinates_class.o
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LDFLAGS) $(LIBS)

clean:
	rm -f *.o core* ruferr
