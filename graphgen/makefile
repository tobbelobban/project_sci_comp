
CC = gcc

GENERATOR_HEADERS = ../graph500-newreference/generator/
GENERATOR_SOURCES = ../graph500-newreference/generator/graph_generator.c ../graph500-newreference/generator/make_graph.c ../graph500-newreference/generator/splittable_mrg.c ../graph500-newreference/generator/utils.c

CFLAGS = -fopenmp -std=gnu11 -g -Wall -I$(GENERATOR_HEADERS)

all: 
	$(CC) $(LDFLAGS) $(CFLAGS) $(GENERATOR_SOURCES) test.c -lm -o a.out