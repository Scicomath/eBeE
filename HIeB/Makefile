testHIeB: testHIeB.o HIeB.o HIeB.hpp HIeB.cpp sq.c sq.h
	g++ -Wall -O -o testHIeB testHIeB.o HIeB.o
testHIeB.o: testHIeB.cpp HIeB.hpp 
	g++ -Wall -O -c testHIeB.cpp
HIeB.o: HIeB.cpp sq.h HIeB.hpp
	g++ -Wall -O -c HIeB.cpp
sq.o: sq.h sq.c 
	g++ -Wall -O -c sq.c
