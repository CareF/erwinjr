cFunctions.so : cFunctions.o
	cc -shared -fPIC -o cFunctions.so cFunctions.o

cFunctions.o : cFunctions.c complex.h
	cc -c -fPIC cFunctions.c

.PHONY : clean
clean :
	rm cFunctions.so cFunctions.o
