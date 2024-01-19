LIBNAME = libkza.so
CFLAGS = -Wall -fPIC
LDFLAGS = -lm
OBJMODULES = kz.o kza.o

.PHONY: clean

all: library

debug: CFLAGS += -DTIMER -DDEBUG
debug: library

kz_threads_strategy_1: CFLAGS += -DKZ_THREADS_STRATEGY_1
kz_threads_strategy_1: LDFLAGS += -lpthread
kz_threads_strategy_1: library

kz_threads_strategy_1_debug: CFLAGS += -DKZ_THREADS_STRATEGY_1 -DTIMER \
									   -DDEBUG
kz_threads_strategy_1_debug: LDFLAGS += -lpthread
kz_threads_strategy_1_debug: library

kz_threads_strategy_2: CFLAGS += -DKZ_THREADS_STRATEGY_1
kz_threads_strategy_2: LDFLAGS += -lpthread
kz_threads_strategy_2: library

kz_threads_strategy_2_debug: CFLAGS += -DKZ_THREADS_STRATEGY_2 -DTIMER \
									   -DDEBUG
kz_threads_strategy_2_debug: LDFLAGS += -lpthread
kz_threads_strategy_2_debug: library

kza_threads_strategy_1: CFLAGS += -DKZA_THREADS_STRATEGY_1
kza_threads_strategy_1: kz_threads_strategy_1

kza_threads_strategy_1_debug: CFLAGS += -DKZA_THREADS_STRATEGY_1 -DTIMER \
									    -DDEBUG
kza_threads_strategy_1_debug: kz_threads_strategy_1_debug

%.o: %.c %.h
	$(CC) $(CFLAGS) -c $< -o $@

library: $(OBJMODULES) 
	$(CC) -shared $^ $(LDFLAGS) -o $(LIBNAME)

clean:
	rm $(LIBNAME) *.o
