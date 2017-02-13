CC = gcc
CXX = g++

CFLAGS = -Wall -O2 -ffast-math
LDLIBS = -larmadillo -lopenblas -llapack

all: main

main: Rcwa.o System.o Mesh.o gauss_legendre.o
	g++ Rcwa.o System.o Mesh.o gauss_legendre.o main.cpp -o main -std=c++11 -O3 -ffast-math -larmadillo -lopenblas -llapack

Mesh.o: Mesh.cpp Rcwa.cpp System.cpp
	g++ -c -std=c++11 Mesh.cpp Rcwa.cpp System.cpp -O3 -ffast-math -larmadillo -lopenblas -llapack

System.o: System.cpp
	g++ -c -std=c++11 System.cpp -O3 -ffast-math -larmadillo -lopenblas -llapack

Rcwa.o: Rcwa.cpp
	g++ -c -std=c++11 Rcwa.cpp -O3 -ffast-math -larmadillo -lopenblas -llapack

gauss_legendre.o: gauss_legendre.c
	gcc -c gauss_legendre.c -O3 -ffast-math -larmadillo -lopenblas -llapack

clean:
	rm *o main
