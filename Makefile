CC=g++ -std=c++11
CFLAGS=-O3 -DNDEBUG

all: cleanall index

test: 
	@$(CC) $(CFLAGS) test.cpp -o heaphull heapchull.a 

HeapCHull.o: HeapCHull.cpp
	$(CC) $(CFLAGS) -c HeapCHull.cpp

index: HeapCHull.o
	ar rc heapchull.a HeapCHull.o

clean:
	-rm *~ *.o *.bak 
cleanall:
	-rm *~ *.o *.bak *.a
