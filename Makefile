LIBNAME = libkza.so
CFLAGS = -Wall -fPIC
LDFLAGS = -lm
OBJMODULES = kz.o kza.o timer.o

.PHONY: clean

all: library

kza_timer: CFLAGS += -DTIMER 
kza_timer: library

prefix_sum: CFLAGS += -DPREFIX_SUM 
prefix_sum: library

prefix_sum-timer: CFLAGS += -DTIMER 
prefix_sum-timer: prefix_sum 

%.o: %.c %.h
	$(CC) $(CFLAGS) -c $< -o $@

library: $(OBJMODULES) 
	$(CC) -shared $^ $(LDFLAGS) -o $(LIBNAME)

clean:
	rm $(LIBNAME) *.o
