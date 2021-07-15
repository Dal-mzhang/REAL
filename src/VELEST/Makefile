FC = gfortran -w
BIN = ../../bin

all: velest clean

velest:velest.o
	$(FC) -o $(BIN)/$@ $^

clean:
	rm -r *.o
