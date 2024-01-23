LIBNAME = libkza.so
CFLAGS = -Wall -fPIC
LDFLAGS = -lm
OBJMODULES = kz.o kza.o

.PHONY: clean

all: library

debug: CFLAGS += -DTIMER -DDEBUG
debug: library

prefix_sum: CFLAGS += -DPREFIX_SUM 
prefix_sum: library

prefix_sum-debug: CFLAGS += -DDEBUG -DTIMER 
prefix_sum-debug: prefix_sum 

threads_client_server: CFLAGS += -DKZ_THREADS_CLIENT_SERVER \
								 -DKZA_THREADS_CLIENT_SERVER
threads_client_server: LDFLAGS += -lpthread
threads_client_server: library

threads_client_server-debug: CFLAGS += -DTIMER -DDEBUG
threads_client_server-debug: threads_client_server

threads_loop: CFLAGS += -DKZ_THREADS_LOOP -DKZA_THREADS_LOOP
threads_loop: LDFLAGS += -lpthread
threads_loop: library

threads_loop-debug: CFLAGS += -DTIMER -DDEBUG
threads_loop-debug: threads_loop

prefix_sum-threads_client_server: CFLAGS += -DPREFIX_SUM
prefix_sum-threads_client_server: threads_client_server

prefix_sum-threads_client_server-debug: CFLAGS += -DPREFIX_SUM
prefix_sum-threads_client_server-debug: threads_client_server-debug

prefix_sum-threads_loop: CFLAGS += -DPREFIX_SUM
prefix_sum-threads_loop: threads_loop

prefix_sum-threads_loop-debug: CFLAGS += -DPREFIX_SUM
prefix_sum-threads_loop-debug: threads_loop-debug

%.o: %.c %.h
	$(CC) $(CFLAGS) -c $< -o $@

library: $(OBJMODULES) 
	$(CC) -shared $^ $(LDFLAGS) -o $(LIBNAME)

clean:
	rm $(LIBNAME) *.o
