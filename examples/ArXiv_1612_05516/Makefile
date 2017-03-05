CFLAGS=-std=c++11 -O3 -march=native
INSTALLPATH=../../
INCLUDES=-I$(INSTALLPATH)src -I$(INSTALLPATH)src/utility -I$(INSTALLPATH)src/mathlib
ARMAINCLUDE=-I$(INSTALLPATH)src/arma -DARMA_DONT_USE_WRAPPER
LIBSPATH=-L$(INSTALLPATH)build -lmesh -L$(ATLASHOME)/lib -L$(LAPACKHOME)/lib -Wl,--start-group -lf77blas -lcblas -latlas -llapack -lcblas -lgfortran -Wl,--end-group
CXX=mpicxx

all:
	$(CXX) main.cpp -o main $(CFLAGS) $(INCLUDES) ${ARMAINCLUDE} $(LIBSPATH) $(LIBS)
