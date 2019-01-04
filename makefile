PROGRAM=pi
CC=mpicc
CFLAG= -Wall -pedantic -std=c99 -I.. -O3 

all:	pi
pi: pi.c
	$(CC) $(CFLAG) -o $(PROGRAM) pi.c
clean:
	rm -f *.o *~ core *.out $(PROGRAM) 
