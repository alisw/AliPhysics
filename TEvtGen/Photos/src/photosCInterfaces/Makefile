include ../make.inc

SOURCES_CC=PhotosParticle.cxx PhotosEvent.cxx PhotosBranch.cxx Photos.cxx
OBJECTS_CC=$(SOURCES_CC:.cxx=.o)

%.o: %.cxx
	$(CC) $(CFLAGS) $(DEBUG) -I. -I../utilities -I../photosFortranInterfaces -I../eventRecordInterfaces -c $< -o $@

all:   $(OBJECTS_CC)

clean:
	rm -f *.o *~;
