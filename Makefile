LIBNAME = libkza.so
CFLAGS = -Wall -fPIC
LDFLAGS = -lm
OBJMODULES = kz.o kza.o

.PHONY: clean

all: library

debug: CFLAGS += -DTIMER -DDEBUG
debug: library

kz_threads_client_server: CFLAGS += -DKZ_THREADS_CLIENT_SERVER
kz_threads_client_server: LDFLAGS += -lpthread
kz_threads_client_server: library

kz_threads_client_server_debug: CFLAGS += -DKZ_THREADS_CLIENT_SERVER -DTIMER \
				-DDEBUG
kz_threads_client_server_debug: LDFLAGS += -lpthread
kz_threads_client_server_debug: library

kz_threads_loop: CFLAGS += -DKZ_THREADS_CLIENT_SERVER
kz_threads_loop: LDFLAGS += -lpthread
kz_threads_loop: library

kz_threads_loop_debug: CFLAGS += -DKZ_THREADS_LOOP -DTIMER \
		       -DDEBUG
kz_threads_loop_debug: LDFLAGS += -lpthread
kz_threads_loop_debug: library

kza_threads_client_server: CFLAGS += -DKZA_THREADS_STRATEGY_1
kza_threads_client_server: kz_threads_client_server

kza_threads_client_server_debug: CFLAGS += -DKZA_THREADS_STRATEGY_1 -DTIMER \
                                 -DDEBUG
kza_threads_client_server_debug: kz_threads_client_server_debug

%.o: %.c %.h
	$(CC) $(CFLAGS) -c $< -o $@

library: $(OBJMODULES) 
	$(CC) -shared $^ $(LDFLAGS) -o $(LIBNAME)

clean:
	rm $(LIBNAME) *.o
