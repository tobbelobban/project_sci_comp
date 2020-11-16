
CC = gcc
CCPP = g++

GENERATOR_HEADERS = ../graph500-newreference/generator/
GENERATOR_SOURCES = ../graph500-newreference/generator/graph_generator.c ../graph500-newreference/generator/make_graph.c ../graph500-newreference/generator/splittable_mrg.c ../graph500-newreference/generator/utils.c

LOCAL_HEADERS = graph.hpp csr.hpp
LOCAL_SOURCES = graph.cpp csr.cpp main.cpp

CPPFLAGS = -fopenmp -std=c++11 -g -Wall -march=native -mavx
CFLAGS = -fopenmp -std=gnu11 -g -Wall -I$(GENERATOR_HEADERS)

all: includes locals
	$(CCPP) -fopenmp -o main.out -I$(GENERATOR_HEADERS) $(wildcard *.o)

includes:
	$(CC) $(CFLAGS) -c $(GENERATOR_SOURCES)

locals:
	$(CCPP) $(CPPFLAGS) -c -I$(GENERATOR_HEADERS) $(LOCAL_SOURCES)

clean:
	rm *.o