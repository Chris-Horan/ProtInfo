all: test

test: test.o
	g++ -o test test.o

test.o: test.cpp Protein.h
	g++ -c test.cpp

clean:
	rm -f test *.o