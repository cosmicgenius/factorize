CC = g++
CFLAGS = -pedantic -Wall -Wextra -lgmp -lgmpxx -g -O3 -pg

siqs: obj/siqs.o obj/sieve_handler.o obj/siever.o obj/util.o obj/gf2.o
	$(CC) -o build/main obj/siqs.o obj/sieve_handler.o obj/siever.o obj/util.o obj/gf2.o $(CFLAGS)

pollard: obj/pollard.o obj/util.o
	$(CC) -o build/main obj/pollard.o obj/util.o $(CFLAGS)

test: obj/test.o obj/util.o obj/gf2.o
	$(CC) -o build/test obj/test.o obj/util.o obj/gf2.o $(CFLAGS)
	./build/test

obj/siqs.o: src/siqs.cpp include/util.hpp include/sieve_handler.hpp include/siever.hpp
	$(CC) -o obj/siqs.o -c src/siqs.cpp $(CFLAGS)

obj/pollard.o: src/pollard.cpp include/util.hpp
	$(CC) -o obj/pollard.o -c src/pollard.cpp $(CFLAGS)

obj/test.o: src/test.cpp include/util.hpp
	$(CC) -o obj/test.o -c src/test.cpp $(CFLAGS)

obj/sieve_handler.o: src/sieve_handler.cpp include/util.hpp include/siever.hpp include/gf2.hpp
	$(CC) -o obj/sieve_handler.o -c src/sieve_handler.cpp $(CFLAGS)

obj/siever.o: src/siever.cpp include/util.hpp
	$(CC) -o obj/siever.o -c src/siever.cpp $(CFLAGS)

obj/util.o: src/util.cpp
	$(CC) -o obj/util.o -c src/util.cpp $(CFLAGS)

obj/gf2.o: src/gf2.cpp
	$(CC) -o obj/gf2.o -c src/gf2.cpp $(CFLAGS)

gen: 
	mkdir -p build
	mkdir -p obj

run: build/main
	./build/main

.PHONY: clean
clean: 
	-rm obj/*.o build/* 

