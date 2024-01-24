CC = g++
CFLAGS = -pedantic -Wall -Wextra -lgmp -lgmpxx -g
OPTFLAGS = -O3 -march=native

siqs: obj/siqs.o obj/sieve_handler.o obj/siever.o obj/util.o obj/gf2.o
	$(CC) -o build/main obj/siqs.o obj/sieve_handler.o obj/siever.o obj/util.o obj/gf2.o $(CFLAGS) $(OPTFLAGS)

pollard: obj/pollard.o obj/util.o
	$(CC) -o build/main obj/pollard.o obj/util.o $(CFLAGS) $(OPTFLAGS)

test: obj/test.o obj/util.o obj/gf2.o
	$(CC) -o build/test obj/test.o obj/util.o obj/gf2.o $(CFLAGS) $(OPTFLAGS)
	./build/test

obj/siqs.o: src/siqs.cpp include/util.hpp include/sieve_handler.hpp include/siever.hpp
	$(CC) -o obj/siqs.o -c src/siqs.cpp $(CFLAGS) $(OPTFLAGS)

obj/pollard.o: src/pollard.cpp include/util.hpp
	$(CC) -o obj/pollard.o -c src/pollard.cpp $(CFLAGS) $(OPTFLAGS)

obj/test.o: src/test.cpp include/util.hpp
	$(CC) -o obj/test.o -c src/test.cpp $(CFLAGS) $(OPTFLAGS)

obj/sieve_handler.o: src/sieve_handler.cpp include/util.hpp include/siever.hpp include/gf2.hpp
	$(CC) -o obj/sieve_handler.o -c src/sieve_handler.cpp $(CFLAGS) $(OPTFLAGS)

obj/siever.o: src/siever.cpp include/util.hpp
	$(CC) -o obj/siever.o -c src/siever.cpp $(CFLAGS) $(OPTFLAGS)

obj/util.o: src/util.cpp
	$(CC) -o obj/util.o -c src/util.cpp $(CFLAGS) $(OPTFLAGS)

obj/gf2.o: src/gf2.cpp
	$(CC) -o obj/gf2.o -c src/gf2.cpp $(CFLAGS) $(OPTFLAGS)

gen: 
	mkdir -p build
	mkdir -p obj

run: build/main
	./build/main

.PHONY: clean
clean: 
	-rm obj/*.o build/* 

