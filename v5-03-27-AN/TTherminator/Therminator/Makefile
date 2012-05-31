CC	 = g++
FF	 = g77
LD	 = g++

BINARIES = therm_events therm_tree
HSOURCES = Parser.cxx ParticleType.cxx DecayChannel.cxx Integrator.cxx ParticleDB.cxx DecayTable.cxx Particle.cxx Event.cxx ParticleDecayer.cxx ReadPar.cxx Hypersurface.cxx
HEADERS  = $(HSOURCES:.cxx=.h)
SOURCES  = $(HSOURCES) therm_events.cxx
OBJECTS  = $(SOURCES:.cxx=.o) 

CPPOPT   = -Wno-deprecated `root-config --cflags`
COPTN    = -c -O0 -g ${CPPOPT} -D_DEBUG_LEVEL_=1
COPT     = -c -O0 -g ${CPPOPT} -D_DEBUG_LEVEL_=1
LLIB     = `root-config --libs` -L$(ROOTSYS)/lib -lm -lgcc -g

all: $(BINARIES)

therm_events: $(OBJECTS)
	$(LD) $(LLIB) -o $@ $(OBJECTS) 

therm_tree: therm_tree.o
	$(LD) $(LLIB) -o $@ therm_tree.o

therm_tree.o: therm_tree.C
	$(CC) $(COPT) -o $@ $^

%.o: %.cxx
	$(CC) $(COPT) -o $@ $< 

%.d: %.cxx
	$(SHELL) -ec '$(CC) -MM $(CPPOPT) $< | sed '\''s/\($*\)\.o[ :]*/\1.o $@ : /g'\'' > $@; [ -s $@ ] || rm -f $@'

include $(SOURCES:.cxx=.d)

clean:
	rm -f *.o therm_events therm_tree therm_events.exe therm_tree.exe *.d 

cleanod:
	rm -f *.o *.d 

package: $(SOURCES) $(HEADERS) Makefile  therminator.in therm_tree.C ./share/particles.data ./share/decays.data THGlobal.h figure*.C runevents.sh
	tar zcvf therminator-2_0.tar.gz $^
