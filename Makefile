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

threads_client_server: CFLAGS += -DKZ_THREADS_CLIENT_SERVER \
								 -DKZA_THREADS_CLIENT_SERVER
threads_client_server: LDFLAGS += -lpthread
threads_client_server: library

threads_client_server-timer: CFLAGS += -DTIMER
threads_client_server-timer: threads_client_server

threads_loop: CFLAGS += -DKZ_THREADS_LOOP -DKZA_THREADS_LOOP
threads_loop: LDFLAGS += -lpthread
threads_loop: library

threads_loop-timer: CFLAGS += -DTIMER
threads_loop-timer: threads_loop

prefix_sum-threads_client_server: CFLAGS += -DPREFIX_SUM
prefix_sum-threads_client_server: threads_client_server

prefix_sum-threads_client_server-timer: CFLAGS += -DPREFIX_SUM
prefix_sum-threads_client_server-timer: threads_client_server-timer

prefix_sum-threads_loop: CFLAGS += -DPREFIX_SUM
prefix_sum-threads_loop: threads_loop

prefix_sum-threads_loop-timer: CFLAGS += -DPREFIX_SUM
prefix_sum-threads_loop-timer: threads_loop-timer

%.o: %.c %.h
	$(CC) $(CFLAGS) -c $< -o $@

library: $(OBJMODULES) 
	$(CC) -shared $^ $(LDFLAGS) -o $(LIBNAME)

clean:
	rm $(LIBNAME) *.o
