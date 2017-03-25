CFLAGS=-std=c++11 -O3
INSTALLPATH=../../
INCLUDES=-I$(INSTALLPATH)src
ARMAINCLUDE=-I$(INSTALLPATH)src/arma -DARMA_DONT_USE_WRAPPER
LIBS=-L$(INSTALLPATH)build -lmesh -lopenblas -llapack
CXX=mpic++

all:
	$(CXX) diagonal.cpp -o diagonal $(CFLAGS) $(INCLUDES) ${ARMAINCLUDE} $(LIBS)
