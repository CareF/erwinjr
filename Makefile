# CC = icc
# CFLAGS = -march=native -Ofast -gcc-name=gcc-6
CC = gcc
CFLAGS = -O2

.PHONY : all
.DEFAULT : all
all: cQCLayers.so cStrata.so cQCLayersMP.so

cQCLayers.so : cQCLayers.o
	$(CC) -shared -fPIC -o cQCLayers.so cQCLayers.o

cQCLayers.o : cQCLayers.c
	$(CC) -c -fPIC $(CFLAGS) cQCLayers.c

cStrata.so : cStrata.o
	$(CC) -shared -fPIC -o cStrata.so cStrata.o

cStrata.o : cStrata.c complex.h
	$(CC) -c -fPIC $(CFLAGS) cStrata.c

cQCLayersMP.so : cQCLayersMP.o
	$(CC) -shared -fPIC -fopenmp -o cQCLayersMP.so cQCLayersMP.o

cQCLayersMP.o : cQCLayers.c
	$(CC) -c -fPIC -fopenmp -D __MP $(CFLAGS) cQCLayers.c -o cQCLayersMP.o

.PHONY : clean
clean :
	rm cQCLayers.so cQCLayers.o cStrata.so cStrata.o cQCLayersMP.so cQCLayersMP.o
