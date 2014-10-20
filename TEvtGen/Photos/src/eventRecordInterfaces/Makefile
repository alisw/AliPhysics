include ../make.inc

SOURCES_CC = PhotosHEPEVTEvent.cxx PhotosHEPEVTParticle.cxx

ifneq ($(HEPMCLOCATION), )
  SOURCES_CC += PhotosHepMCParticle.cxx PhotosHepMCEvent.cxx
endif

OBJECTS_CC=$(SOURCES_CC:.cxx=.o)

%.o: %.cxx
	$(CC) $(CFLAGS) $(DEBUG) -I$(HEPMCLOCATION)/include -I. -I../utilities -I../photosFortranInterfaces -I../photosCInterfaces -c $< -o $@

all:   $(OBJECTS_CC)

clean:
	rm -f *.o *~;
