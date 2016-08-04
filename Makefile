CC=gcc
CFLAGS=-O3 -std=c99

objects = eBeE.o eBFun.o rhoFun.o sq.o sqrtStoY.o mean_eB.o

eBeE: $(objects)
	$(CC) $(CFLAGS) -o eBeE $(objects) -I. libcuba.a -lgsl -lgslcblas -lm

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

mean_eB.o: mean_eB.c udStruct.h eBFun.h sq.h
	$(CC) $(CFLAGS) -c mean_eB.c

.PHONY: clean
clean:
	-rm eBeE $(objects)
