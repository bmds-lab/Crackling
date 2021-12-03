# the compiler:
CC = g++

# compiler flags:
CFLAGS = -O3 -std=c++11 -fopenmp -mpopcnt

# define any directories containing header files other than /usr/include
INCLUDES = -Iparallel_hashmap

all : isslScoreOfftargets isslCreateIndex

isslScoreOfftargets : isslScoreOfftargets.cpp
	$(CC) $(CFLAGS) $(INCLUDES) -o $@ $^

isslCreateIndex : isslCreateIndex.cpp
	$(CC) $(CFLAGS) $(INCLUDES) -o $@ $^

clean:
	$(RM) isslScoreOfftargets isslCreateIndex
