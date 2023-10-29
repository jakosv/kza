LIBNAME = libkza.so
CFLAGS = -Wall -fPIC
LDFLAGS = -lm
OBJMODULES = kz.o kza.o

.PHONY: clean

all: library

debug: CFLAGS += -DTIMER
debug: library

kz_parallel: CFLAGS += -DKZ_PARALLEL
kz_parallel: LDFLAGS += -lpthread
kz_parallel: library

kz_parallel_debug: CFLAGS += -DKZ_PARALLEL -DTIMER
kz_parallel_debug: LDFLAGS += -lpthread
kz_parallel_debug: library

speed_hacks: CFLAGS += -DSPEED library

%.o: %.c %.h
	$(CC) $(CFLAGS) -c $< -o $@

library: $(OBJMODULES) 
	$(CC) -shared $^ $(LDFLAGS) -o $(LIBNAME)

clean:
	rm $(LIBNAME) *.o
