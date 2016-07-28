CC=gcc
CFLAGS=-O3 -std=c99

objects = eBeE.o eBFun.o rhoFun.o sq.o sqrtStoY.o

eBeE: $(objects)
	$(CC) $(CFLAGS) -o eBeE $(objects) -I. libcuba.a -lgsl -lgslcblas -lm
testrho: testrho.o rhoFun.o sq.o
	$(CC) $(CFLAGS) -o testrho testrho.o rhoFun.o sq.o -I. libcuba.a -lgsl -lgslcblas -lm

testrho.o: testrho.c udStruct.h rhoFun.h sq.h
	$(CC) $(CFLAGS) -c testrho.c

eBeE.o: eBeE.c udStruct.h eBFun.h sqrtStoY.h
	$(CC) $(CFLAGS) -c eBeE.c

eBFun.o: eBFun.c cuba.h udStruct.h rhoFun.h eBFun.h sq.h
	$(CC) $(CFLAGS) -c eBFun.c

rhoFun.o: rhoFun.c rhoFun.h sq.h
	$(CC) $(CFLAGS) -c rhoFun.c

sq.o: sq.c sq.h
	$(CC) $(CFLAGS) -c sq.c

sqrtStoY.o: sqrtStoY.c sqrtStoY.h sq.h
	$(CC) $(CFLAGS) -c sqrtStoY.c

.PHONY: clean
clean:
	-rm eBeE testrho $(objects)
