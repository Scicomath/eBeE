eBeE: eBeE.o eBFun.o rhoFun.o sq.o
	gcc -O3 -o eBeE eBeE.o eBFun.o rhoFun.o sq.o -I. libcuba.a -lgsl -lgslcblas -lm

eBeE.o: eBeE.c udStruct.h eBFun.h
	gcc -O3 -c eBeE.c

eBFun.o: eBFun.c cuba.h udStruct.h rhoFun.h eBFun.h sq.h
	gcc -O3 -c eBFun.c

rhoFun.o: rhoFun.c rhoFun.h sq.h
	gcc -O3 -c rhoFun.c

sq.o: sq.c sq.h
	gcc -O3 -c sq.c
