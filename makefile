CC=mpicc
CXX=mpicxx
LINK=mpicxx
CFLAGS=-I.

%.o: %.c
				$(CC) -c -o $@ $< $(CFLAGS)

%.o: %.cpp
				$(CXX) -c -o $@ $< $(CFLAGS)

all: samplesort seqsort

seqsort: seqsort.o 
			  $(LINK) -o $@ $< 

samplesort: samplesort.o 
			  $(LINK) -o $@ $< 
#data: data.o 
#			  $(LINK) -o $@ $< 


#select: select.o 
#			  $(LINK) -o $@ $< 

clean:
	rm -f *.o *~ seqsort samplesort