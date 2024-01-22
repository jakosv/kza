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

kz_threads_client_server_debug: CFLAGS += -DTIMER -DDEBUG
kz_threads_client_server_debug: kz_threads_client_server

kz_threads_loop: CFLAGS += -DKZ_THREADS_LOOP
kz_threads_loop: LDFLAGS += -lpthread
kz_threads_loop: library

kz_threads_loop_debug: CFLAGS += -DTIMER -DDEBUG
kz_threads_loop_debug: kz_threads_loop

kza_threads_client_server: CFLAGS += -DKZA_THREADS_CLIENT_SERVER
kza_threads_client_server: kz_threads_client_server

kza_threads_client_server_debug: CFLAGS += -DTIMER -DDEBUG
kza_threads_client_server_debug: kz_threads_client_server_debug

kza_threads_loop: CFLAGS += -DKZA_THREADS_LOOP
kza_threads_loop: kz_threads_loop

kza_threads_loop_debug: CFLAGS += -DKZA_THREADS_LOOP -DTIMER -DDEBUG
kza_threads_loop_debug: kz_threads_loop_debug

%.o: %.c %.h
	$(CC) $(CFLAGS) -c $< -o $@

library: $(OBJMODULES) 
	$(CC) -shared $^ $(LDFLAGS) -o $(LIBNAME)

clean:
	rm $(LIBNAME) *.o
