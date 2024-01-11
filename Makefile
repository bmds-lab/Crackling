# the compiler:
CC = g++

# compiler flags:
CFLAGS = -O3 -std=c++11 -fopenmp -mpopcnt

# define any directories containing header files other than /usr/include
INCLUDES = -Isrc/ISSL/include

all : bin_check isslScoreOfftargets isslCreateIndex

bin_check:
	@if [ ! -d bin ]; then \
		mkdir bin; \
	fi

isslScoreOfftargets : src/ISSL/isslScoreOfftargets.cpp
	$(CC) $(CFLAGS) $(INCLUDES) -o bin/$@ $^

isslCreateIndex : src/ISSL/isslCreateIndex.cpp
	$(CC) $(CFLAGS) $(INCLUDES) -o bin/$@ $^

clean:
	$(RM) bin/isslScoreOfftargets bin/isslCreateIndex
