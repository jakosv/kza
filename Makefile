LIBNAME = libkza
CFLAGS = -Wall -fPIC #-DPARALLEL -DTIMER -DSPEED
LDFLAGS = -lm -lpthread
OBJMODULES = kz.o kza.o

.PHONY: clean

%.o: %.c %.h
	$(CC) $(CFLAGS) -c $< -o $@

$(LIBNAME).so: $(OBJMODULES) 
	$(CC) -shared $(CFALGS) $^ $(LDFLAGS) -o $(LIBNAME).so

clean:
	rm $(LIBNAME).so *.o
