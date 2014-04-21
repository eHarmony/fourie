
LDFLAGS := $(LDFLAGS) -L/usr/local/lib

LDFLAGS := $(LDFLAGS) -lm -l boost_program_options -l z -l pthread
CFLAGS := $(CFLAGS) -Wall -pedantic

BINARIES := fourie

all: fourie

fourie: fourie.o
	$(CXX) fourie.o $(LDFLAGS) -o fourie

fourie.o: fourie.cpp
	$(CXX)  $(CFLAGS) -c fourie.cpp -o fourie.o

clean:
	rm -f  *.o $(BINARIES) *~

