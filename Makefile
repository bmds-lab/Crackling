# the compiler:
CC = g++

# compiler flags:
CFLAGS = -O3 -std=c++11 -fopenmp -mpopcnt

# define any directories containing header files other than /usr/include
INCLUDES = -Iparallel_hashmap

# define the source files
SRC = isslScoreOfftargets.cpp isslCreateIndex.cpp

# the build target executables:
TARGET = isslScoreOfftargets isslCreateIndex

all:	$(TARGET)

$(TARGET) : $(SRC) 
	$(CC) $(CFLAGS) $(INCLUDES) -o $@ $<

clean:
	$(RM) $(TARGET)

depend: $(SRC)
	makedepend $(INCLUDES) $^

# DO NOT DELETE THIS LINE -- make depend needs it