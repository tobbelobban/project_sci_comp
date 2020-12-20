
CC = gcc
CCPP = g++

LOCAL_HEADERS = *.hpp
LOCAL_SOURCES = *.cpp

CPPFLAGS = -fopenmp -std=c++11 -g -Wall -march=native -mavx

all: $(LOCAL_HEADERS) $(LOCAL_SOURCES)
	$(CCPP) -o bfs.out $(CPPFLAGS) $(LOCAL_SOURCES)

clean:
	rm *.o