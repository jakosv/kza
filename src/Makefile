LIBNAME = libkza.so
CXXFLAGS = -std=c++20 -fPIC -O2 -Wall
LDFLAGS = -lstdc++
OBJMODULES = kz.o kza.o

.PHONY: clean

all: library

debug: CFLAGS += -DDEBUG 
debug: CXXFLAGS += -DDEBUG 
debug: library

kza_timer: CXXFLAGS += -DTIMER 
kza_timer: library

prefix_sum: CXXFLAGS += -DPREFIX_SUM 
prefix_sum: library

prefix_sum-timer: CXXFLAGS += -DTIMER 
prefix_sum-timer: prefix_sum 

%.o: %.cpp %.h %.hpp
	$(CC) $(CXXFLAGS) -c $< -o $@

library: $(OBJMODULES) 
	$(CC) -shared $^ $(LDFLAGS) -o $(LIBNAME)

clean:
	rm $(LIBNAME) *.o
