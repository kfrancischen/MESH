CFLAGS=-std=c++11 -O3
INCLUDES=-I../../src -I../../src/utility -I../../src/mathlib
LIBS=-L../../build -lmesh -larmadillo -lopenblas -llapack
CXX=mpic++

all:
	$(CXX) $(CFLAGS) $(INCLUDES) test.cpp -o test $(LIBS)
