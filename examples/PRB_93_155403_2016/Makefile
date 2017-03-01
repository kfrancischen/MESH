CFLAGS=-std=c++11 -O3
INSTALLPATH=../../
INCLUDES=-I$(INSTALLPATH)src -I$(INSTALLPATH)src/utility -I$(INSTALLPATH)src/mathlib
ARMAINCLUDE=-I$(INSTALLPATH)src/arma -DARMA_DONT_USE_WRAPPER
LIBS=-L$(INSTALLPATH)build -lmesh -lopenblas -llapack
CXX=mpic++

all:
	$(CXX) $(CFLAGS) $(INCLUDES) ${ARMAINCLUDE} main.cpp -o main $(LIBS)
