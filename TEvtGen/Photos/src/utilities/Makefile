include ../make.inc

SOURCES_CC=Log.cxx PhotosRandom.cxx PhotosDebugRandom.cxx PhotosUtilities.cxx
OBJECTS_CC=$(SOURCES_CC:.cxx=.o)

%.o: %.cxx
	$(CC) $(CFLAGS) $(DEBUG) -I. -I../photosCInterfaces -I../photosFortranInterfaces -I../eventRecordInterfaces -c $< -o $@

all:   $(OBJECTS_CC)

clean:
	rm -f *.o *~;
