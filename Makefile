CC = g++
CFLAG = -std=c++17 -pthread
SRC = path.cpp

all: $(SRC)
	$(CC) $(CFLAG) $(SRC) -o lateracer

.PHONY: clean
clean:
	rm -rf lateracer
