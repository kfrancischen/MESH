CFLAGS=-std=c++11 -O3 -ffast-math -march=native -fopenmp
INSTALLPATH=../../
INCLUDES=-I$(INSTALLPATH)/src
ARMAINCLUDE=-I$(INSTALLPATH)/src/arma -DARMA_DONT_USE_WRAPPER
LIBS=-L$(INSTALLPATH)build -lmesh -lopenblas -llapack
CXX=g++

all:
	$(CXX) $(CFLAGS) $(INCLUDES) ${ARMAINCLUDE} test.cpp -o test $(LIBS)
