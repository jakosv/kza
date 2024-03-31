LIBNAME = libkza.so
CXXFLAGS = -std=c++20 -Wall -fPIC -O2
LDFLAGS = -lstdc++
OBJMODULES = kz.o kza.o timer.o

.PHONY: clean

all: library

debug: CFLAGS += -DDEBUG 
debug: CXXFLAGS += -DDEBUG 
debug: library

kza_timer: CFLAGS += -DTIMER 
kza_timer: CXXFLAGS += -DTIMER 
kza_timer: library

prefix_sum: CFLAGS += -DPREFIX_SUM 
prefix_sum: CXXFLAGS += -DPREFIX_SUM 
prefix_sum: library

prefix_sum-timer: CFLAGS += -DTIMER 
prefix_sum-timer: CXXFLAGS += -DTIMER 
prefix_sum-timer: prefix_sum 

%.o: %.cpp %.h %.hpp
	$(CC) $(CXXFLAGS) -c $< -o $@

library: $(OBJMODULES) 
	$(CC) -shared $^ $(LDFLAGS) -o $(LIBNAME)

clean:
	rm $(LIBNAME) *.o
