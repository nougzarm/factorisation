PROG = test

LIBR = factorisation

OBJS = test.o $(LIBR).o
CC = gcc
CPFLAGS = -Wall -Wpointer-arith -O2


test: test.o $(LIBR).o
	$(CC) $(CPFLAGS) -o test  test.o $(LIBR).o -lm

test.o: $(LIBR).h test.h test.c
	$(CC) $(CPFLAGS) -c test.c -L/usr/local/lib/ -lgmp

$(LIBR).o: $(LIBR).c $(LIBR).h
	$(CC) $(CPFLAGS) -c $(LIBR).c -L/usr/local/lib/ -lgmp

clean:
	rm -f $(OBJS)