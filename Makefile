cFunctions.so : cFunctions.o
	cc -shared -fPIC -o cFunctions.so cFunctions.o

cFunctions.o : cFunctions.c complex.h
	cc -c -fPIC cFunctions.c

cFunctionsMP.so : cFunctionsMP.o
	cc -shared -fPIC -fopenmp -o cFunctionsMP.so cFunctionsMP.o

cFunctionsMP.o : cFunctions.c complex.h
	cc -c -fPIC -fopenmp -D __MP cFunctions.c -o cFunctionsMP.o

.PHONY : clean
clean :
	rm cFunctions.so cFunctions.o cFunctionsMP.so cFunctionsMP.o
