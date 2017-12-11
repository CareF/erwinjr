# CC = icc
# CFLAGS = -march=native -Ofast -gcc-name=gcc-6 -Wall -c -fPIC
CC = gcc
CFLAGS = -O2 -Wall -c -fPIC

.PHONY : all
.DEFAULT : all
all: cQCLayers.so cStrata.so cQCLayersMP.so

cQCLayers.so : cQCLayers.o
	$(CC) -shared -fPIC $< -o $@ 

cQCLayers.o : cQCLayers.c
	$(CC) $(CFLAGS) $<

cStrata.so : cStrata.o
	$(CC) -shared -fPIC $< -o $@ 

cStrata.o : cStrata.c complex.h
	$(CC) $(CFLAGS) $<

cQCLayersMP.so : cQCLayersMP.o
	$(CC) -shared -fPIC -fopenmp $< -o $@ 

cQCLayersMP.o : cQCLayers.c
	$(CC) -fopenmp -D __MP $(CFLAGS) $< -o $@

.PHONY : clean
clean :
	rm cQCLayers.so cQCLayers.o cStrata.so cStrata.o cQCLayersMP.so cQCLayersMP.o
