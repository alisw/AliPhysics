include make.inc

ifneq ($(PHOTOSLOCATION), )
include $(PHOTOSLOCATION)/make.inc
endif

ROOTLIB       :=  $(shell root-config --glibs)
ROOTINC       :=  $(shell root-config --cflags)
ROOTLIBPATH   :=  $(shell root-config --libdir)
TAUOLA_LIBS    = -L$(TAUOLALOCATION)/lib -lTauolaCxxInterface -lTauolaFortran
PYTHIA_LIBS    = -L$(PYTHIALOCATION)/lib/archive -lpythia8 -llhapdfdummy -lhepmcinterface
MC_TESTER_LIBS = -L$(MCTESTERLOCATION)/lib -lHEPEvent -lHepMCEvent -lMCTester
PHOTOS_LIBS    = -L$(PHOTOSLOCATION)/lib -lPhotosCxxInterface -lPhotosFortran $(HEPMC_LIBS)
HEPMC_LIBS     =

ifneq ($(HEPMCLOCATION), )
  override HEPMC_LIBS += -L$(HEPMCLOCATION)/lib -lHepMC
endif

all:
	@make photos_hepevt_example.exe
ifneq ($(HEPMCLOCATION), )
	@make photos_standalone_example.exe
ifneq ($(PYTHIALOCATION), )
	@make single_photos_gun_example.exe
	@make photosLCG_pythia_example.exe
ifneq ($(MCTESTERLOCATION), )
	@make photos_pythia_example.exe
ifneq ($(TAUOLALOCATION), )
	@make tauola_photos_pythia_example.exe
endif
endif
endif
endif
	@echo "##################################################################"
	@echo " Examples compiled successfully. Available examples:"
	@if test -e photos_hepevt_example.exe; then \
	echo " ./photos_hepevt_example.exe"; fi
	@if test -e photos_standalone_example.exe; then \
	echo " ./photos_standalone_example.exe"; fi
	@if test -e single_photos_gun_example.exe; then \
	echo " ./single_photos_gun_example.exe"; fi
	@if test -e photosLCG_pythia_example.exe; then \
	echo " ./photosLCG_pythia_example.exe"; fi
	@if test -e photos_pythia_example.exe; then \
	echo " ./photos_pythia_example.exe"; fi
	@if test -e tauola_photos_pythia_example.exe; then \
	echo " ./tauola_photos_pythia_example.exe"; fi
	@echo "##################################################################"

# Universal example (all paths presents)
%.o: %.cxx
	$(CC) $(CFLAGS) -I$(PHOTOSLOCATION)/include -I$(TAUOLALOCATION)/include -I$(PYTHIALOCATION)/include \
	-I$(MCTESTERLOCATION)/include -I$(HEPMCLOCATION)/include $(ROOTINC) -c $<

# Universal example (all paths presents) - linking procedure
%.exe: %.o
	$(CC) -o $@ $(LDFLAGS) $< $(ROOTLIB) \
	$(TAUOLA_LIBS) $(PYTHIA_LIBS) $(MC_TESTER_LIBS) $(HEPMC_LIBS) $(PHOTOS_LIBS) \
	-Wl,-rpath,$(PREFIX)/lib \
	-Wl,-rpath,$(PHOTOSLOCATION)/lib \
	-Wl,-rpath,$(ROOTLIBPATH) \
	-Wl,-rpath,$(TAUOLALOCATION)/lib \
	-Wl,-rpath,$(HEPMCLOCATION)/lib \
	-Wl,-rpath,$(PYTHIALOCATION)/lib \
	-Wl,-rpath,$(MCTESTERLOCATION)/lib

# Version for pythia example - HepMC, Pythia8 and MC-Tester required
photos_pythia_example.exe: photos_pythia_example.cxx
	$(CC) $(CFLAGS) -I$(PHOTOSLOCATION)/include -I$(PYTHIALOCATION)/include \
	-I$(MCTESTERLOCATION)/include -I$(HEPMCLOCATION)/include $(ROOTINC) -c $<
	$(CC) -o $@ $(LDFLAGS) $(<:.cxx=.o) $(ROOTLIB) \
	$(PYTHIA_LIBS) $(MC_TESTER_LIBS) $(HEPMC_LIBS) $(PHOTOS_LIBS) \
	-Wl,-rpath,$(PREFIX)/lib \
	-Wl,-rpath,$(PHOTOSLOCATION)/lib \
	-Wl,-rpath,$(ROOTLIBPATH) \
	-Wl,-rpath,$(HEPMCLOCATION)/lib \
	-Wl,-rpath,$(PYTHIALOCATION)/lib \
	-Wl,-rpath,$(MCTESTERLOCATION)/lib

# Version for photos gun example - HepMC and Pythia8 required
single_photos_gun_example.exe: single_photos_gun_example.cxx
	$(CC) $(CFLAGS) -I$(PHOTOSLOCATION)/include -I$(PYTHIALOCATION)/include \
	-I$(HEPMCLOCATION)/include -c $<
	$(CC) -o $@ $(LDFLAGS) $(<:.cxx=.o) $(HEPMC_LIBS) $(PYTHIA_LIBS) $(PHOTOS_LIBS) \
	-Wl,-rpath,$(PREFIX)/lib \
	-Wl,-rpath,$(PHOTOSLOCATION)/lib \
	-Wl,-rpath,$(HEPMCLOCATION)/lib \
	-Wl,-rpath,$(PYTHIALOCATION)/lib

# Version for photos LCG example - HepMC and Pythia8 required
photosLCG_pythia_example.exe: photosLCG_pythia_example.cxx
	$(CC) $(CFLAGS) -I$(PHOTOSLOCATION)/include -I$(PYTHIALOCATION)/include \
	-I$(HEPMCLOCATION)/include -c $<
	$(CC) -o $@ $(LDFLAGS) $(<:.cxx=.o) $(HEPMC_LIBS) $(PYTHIA_LIBS) $(PHOTOS_LIBS) \
	-Wl,-rpath,$(PREFIX)/lib \
	-Wl,-rpath,$(PHOTOSLOCATION)/lib \
	-Wl,-rpath,$(HEPMCLOCATION)/lib \
	-Wl,-rpath,$(PYTHIALOCATION)/lib

# Version for photos standalone example - only HepMC required
photos_standalone_example.exe: photos_standalone_example.cxx
	$(CC) $(CFLAGS) -I$(PHOTOSLOCATION)/include -I$(HEPMCLOCATION)/include -c $<
	$(CC) -o $@ $(LDFLAGS) $(<:.cxx=.o) $(HEPMC_LIBS) $(PHOTOS_LIBS) \
	-Wl,-rpath,$(PREFIX)/lib \
	-Wl,-rpath,$(PHOTOSLOCATION)/lib \
	-Wl,-rpath,$(HEPMCLOCATION)/lib

# Version for photos HEPEVT example - no other libraries are required
photos_hepevt_example.exe: photos_hepevt_example.f photos_hepevt_example_interface.cxx
	$(F77) $(FFLAGS) -c $<
	$(CC)  $(CFLAGS) -I$(PHOTOSLOCATION)/include -c photos_hepevt_example_interface.cxx
	$(F77) -o $@ $(LDFLAGS) $(<:.f=.o) photos_hepevt_example_interface.o $(HEPMC_LIBS) $(PHOTOS_LIBS) \
	-Wl,-rpath,$(PREFIX)/lib \
	-Wl,-rpath,$(PHOTOSLOCATION)/lib

clean:
	rm -f *.o *.a *~ *.exe

Clean: clean
	rm -rf config.log config.status autom4te.cache configure.paths.sh
	rm -f make.inc

make.inc:
	@echo ""
	@echo "Please execute ./configure first!"
	@echo ""
	@false

$(PHOTOSLOCATION)/make.inc:
	@echo ""
	@echo "Please compile the interface first!"
	@echo ""
	@false
