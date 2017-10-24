# CC = icc
# CFLAGS = -march=native -O3 -gcc-name=gcc-6
CC = gcc
CFLAGS = -march=native -Ofast
cFunctions.so : cFunctions.o
	$(CC) -shared -fPIC -o cFunctions.so cFunctions.o

cFunctions.o : cFunctions.c complex.h
	$(CC) -c -fPIC $(CFLAGS) cFunctions.c

cFunctionsMP.so : cFunctionsMP.o
	$(CC) -shared -fPIC -fopenmp -o cFunctionsMP.so cFunctionsMP.o

cFunctionsMP.o : cFunctions.c complex.h
	$(CC) -c -fPIC -fopenmp -D __MP $(CFLAGS) cFunctions.c -o cFunctionsMP.o

.PHONY : clean
clean :
	rm cFunctions.so cFunctions.o cFunctionsMP.so cFunctionsMP.o
