CC = g++
CFLAGS = -pedantic -Wall -Wextra -lgmp -g -O3

pollard: obj/pollard.o obj/util.o
	$(CC) -o build/main obj/pollard.o obj/util.o $(CFLAGS)

obj/pollard.o: src/pollard.cpp include/util.hpp
	$(CC) -o obj/pollard.o -c src/pollard.cpp $(CFLAGS)

obj/util.o: src/util.cpp
	$(CC) -o obj/util.o -c src/util.cpp $(CFLAGS)

gen: 
	mkdir -p build
	mkdir -p obj

run: build/main
	@./build/main

.PHONY: clean
clean: 
	-rm obj/*.o build/*

