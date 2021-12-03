# the compiler:
CC = g++

# compiler flags:
CFLAGS = -O3 -std=c++11 -fopenmp -mpopcnt

# define any directories containing header files other than /usr/include
INCLUDES = -Icpp/parallel_hashmap

all : isslScoreOfftargets isslCreateIndex

isslScoreOfftargets : ./cpp/isslScoreOfftargets.cpp
	$(CC) $(CFLAGS) $(INCLUDES) -o ./bin/$@ $^

isslCreateIndex : ./cpp/isslCreateIndex.cpp
	$(CC) $(CFLAGS) $(INCLUDES) -o ./bin/$@ $^

clean:
	$(RM) isslScoreOfftargets isslCreateIndex
