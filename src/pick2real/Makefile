CC = gcc -Os -mcmodel=medium -w
LIBS = -lm
BIN= ../../bin

pick2real: pick2real.o
	$(CC) -o $(BIN)/$@ $^ $(LIBS)

clean:
	rm -f *.o
