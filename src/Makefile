all: main

main: Rcwa.o
	g++ Rcwa.o main.cpp -o main -O2 -larmadillo

main.o: main.cpp
	g++ -c main.cpp -O2 -larmadillo

Rcwa.o: Rcwa.cpp
	g++ -c Rcwa.cpp -O2 -larmadillo

clean:
	rm *o main
