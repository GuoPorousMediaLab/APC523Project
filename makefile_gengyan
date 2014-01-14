CC = g++
BIN = project
OBJ = main.o model.o cell.o interface.o slope_limiter.o riemann_solver.o io.o grid.o
DEBUG = -g
PROF = -pg
CFLAG = -Wall -c $(DEBUG)
LFLAG = -Wall $(DEBUG)

all: $(BIN)

$(BIN): $(OBJ)
	$(CC) $(LFLAG) -o $@ $^

%.o: %.cpp
	$(CC) $(CFLAG) -o $@ $<

clean:
	rm -f $(OBJ) $(BIN) *.txt *~
