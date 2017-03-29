CFLAGS=-std=c++11 -O3 -ffast-math -march=native -fopenmp
MESHPATH=../../
INCLUDES=-I$(MESHPATH)/src
ARMAINCLUDE=-I$(MESHPATH)/src/arma -DARMA_DONT_USE_WRAPPER
LIBS=-L$(MESHPATH)/build -lmesh -lopenblas -llapack
CXX=g++

all:
	$(CXX) $(CFLAGS) $(INCLUDES) ${ARMAINCLUDE} test.cpp -o main $(LIBS)
