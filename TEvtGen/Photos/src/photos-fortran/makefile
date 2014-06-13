include ../make.inc

SOURCES_C=photosC.cxx forZ-MEc.cxx forW-MEc.cxx
OBJECTS_C=$(SOURCES_C:.cxx=.o)

%.o: %.cxx
	$(CC) $(CFLAGS) $(DEBUG) -I. -I../photosCInterfaces -I../photosFortranInterfaces -I../utilities -c $< -o $@

all: $(OBJECTS_C)

clean:
	rm -f *.o *~;
