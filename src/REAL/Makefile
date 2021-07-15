#MAC openMP
CC = gcc-10 -Os -fopenmp -W 

#Linux openMP
#CC = gcc -Os  -fopenmp -W 

#Serial
#CC = gcc -Os -W 

LIBS = -lm 

all: REAL clean

BIN = ../../bin

REAL:REAL.o
	$(CC) -o $(BIN)/$@ $^ $(LIBS)

clean:
	rm -r *.o
