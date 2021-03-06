# CC = icc
# CFLAGS = -march=native -Ofast -gcc-name=gcc-6 -Wall -c -fPIC
CC = gcc
CFLAGS = -O2 -Wall -fPIC

.PHONY : all
.DEFAULT : all
all: cQCLayers.so cStrata.so 

cQCLayers.so : cQCLayers.o
	$(CC) -shared -fPIC $< -o $@ 

cQCLayers.o : cQCLayers.c
	$(CC) $(CFLAGS) -c $<

cStrata.so : cStrata.o
	$(CC) -shared -fPIC $< -o $@ 

cStrata.o : cStrata.c complex.h
	$(CC) $(CFLAGS) -c $<

cQCLayersMP.so : cQCLayersMP.o
	$(CC) -shared -fPIC -fopenmp $< -o $@ 

cQCLayersMP.o : cQCLayers.c
	$(CC) -fopenmp -D __MP $(CFLAGS) -c $< -o $@

.PHONY : clean
clean :
	rm cQCLayers.so cQCLayers.o cStrata.so cStrata.o cQCLayersMP.so cQCLayersMP.o
