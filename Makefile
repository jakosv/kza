LIBNAME = libkza
CFLAGS = -Wall -fPIC
LDFLAGS = -lm
OBJMODULES = kz.o kza.o

.PHONY: clean

%.o: %.c %.h
	$(CC) $(CFLAGS) -c $< -o $@

$(LIBNAME).so: $(OBJMODULES) 
	$(CC) -shared $(CFALGS) $^ $(LDFLAGS) -o $(LIBNAME).so

clean:
	rm $(LIBNAME).so *.o
