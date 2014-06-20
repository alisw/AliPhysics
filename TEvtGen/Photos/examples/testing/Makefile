include ../make.inc
include $(PHOTOSLOCATION)/make.inc

ROOTLIB      := $(shell root-config --glibs)
ROOTINC      := $(shell root-config --cflags)
ROOTLIBPATH  := $(shell root-config --libdir)
PYTHIA_LIBS    = -L$(PYTHIALOCATION)/lib/archive -lpythia8 -llhapdfdummy -lhepmcinterface
MC_TESTER_LIBS = -L$(MCTESTERLOCATION)/lib -lHEPEvent -lHepMCEvent -lMCTester
HEPMC_LIBS     = -L$(HEPMCLOCATION)/lib -lHepMC
TAUOLA_LIBS    = -L$(TAUOLALOCATION)/lib -lTauolaCxxInterface -lTauolaFortran
PHOTOS_LIBS    = -L$(PHOTOSLOCATION)/lib -lPhotosCxxInterface -lPhotosFortran $(HEPMC_LIBS)

all:
ifneq ($(HEPMCLOCATION), )
ifneq ($(PYTHIALOCATION), )
ifneq ($(MCTESTERLOCATION), )
	@make photos_test.exe
ifneq ($(TAUOLALOCATION), )
	@make photos_tauola_test.exe
endif
endif
endif
endif
	@echo "##################################################################"
	@echo " Available programs:"
	@if test -e photos_test.exe; then \
	  echo " ./photos_test.exe"; else \
	  echo ""; \
	  echo " WARNING: Main program used for tests - 'photos_test.exe' "; \
	  echo "          requires HepMC, Pythia8 and MC-TESTER."; \
	fi
	@if test -e photos_tauola_test.exe; then \
	  echo " ./photos_tauola_test.exe"; else \
	  echo ""; \
	  echo " WARNING: Second program used for tests - 'photos_tauola_test.exe' "; \
	  echo "          requires HepMC, Tauola++, Pythia8 and MC-TESTER."; \
	fi
	@echo ""
	@echo "   To execute all available tests: 'make run',"
	@echo "   For single test: cd to chosen sub-directory and 'make' "
	@echo "##################################################################"

%.exe: %.cxx
	$(CC) $(CFLAGS) -I$(PHOTOSLOCATION)/include -I$(HEPMCLOCATION)/include \
        $(LDFLAGS) $(HEPMC_LIBS) $(PHOTOS_LIBS) $< -o $@ \
        -Wl,-rpath,$(PHOTOSLOCATION)/lib \
        -Wl,-rpath,$(HEPMCLOCATION)/lib

photos_test.exe: photos_test.cxx
	$(CC) $(CFLAGS) -I$(PHOTOSLOCATION)/include -I$(PYTHIALOCATION)/include -I$(MCTESTERLOCATION)/include -I$(HEPMCLOCATION)/include \
	$(ROOTINC) $< -o $@ \
	$(LDFLAGS) $(ROOTLIB) $(PYTHIA_LIBS) $(MC_TESTER_LIBS) $(HEPMC_LIBS) $(PHOTOS_LIBS) \
	-Wl,-rpath,$(PHOTOSLOCATION)/lib \
	-Wl,-rpath,$(ROOTLIBPATH) \
	-Wl,-rpath,$(HEPMCLOCATION)/lib \
	-Wl,-rpath,$(PYTHIALOCATION)/lib \
	-Wl,-rpath,$(MCTESTERLOCATION)/lib

photos_tauola_test.exe: photos_tauola_test.cxx
	$(CC) $(CFLAGS) -I$(PHOTOSLOCATION)/include -I$(TAUOLALOCATION)/include -I$(PYTHIALOCATION)/include -I$(MCTESTERLOCATION)/include \
	-I$(HEPMCLOCATION)/include $(ROOTINC) $< -o $@ \
	$(LDFLAGS) $(ROOTLIB) $(PYTHIA_LIBS) $(MC_TESTER_LIBS) $(HEPMC_LIBS) $(TAUOLA_LIBS) $(PHOTOS_LIBS) \
	-Wl,-rpath,$(PHOTOSLOCATION)/lib \
	-Wl,-rpath,$(TAUOLALOCATION)/lib \
	-Wl,-rpath,$(ROOTLIBPATH) \
	-Wl,-rpath,$(HEPMCLOCATION)/lib \
	-Wl,-rpath,$(PYTHIALOCATION)/lib \
	-Wl,-rpath,$(MCTESTERLOCATION)/lib

make.inc:
	@echo ""
	@echo "Please execute ./configure in main examples directory first!"
	@echo ""
	@false

$(PHOTOSLOCATION)/make.inc:
	@echo ""
	@echo "Please compile the interface first!"
	@echo ""
	@false

run:
	make -C Zee
	make -C Zmumu
	make -C ZmumuNLO
	make -C Wenu
	make -C Wmunu
	make -C WmunuNLO
	make -C ttbar
ifneq ($(TAUOLALOCATION), )
	make -C Htautau
	make -C Ztautau
	make -C ScalNLO
endif

clean:
	rm -f *.o *~ *.exe
	rm -f *~ */*~

clobber: clean
	make -C Zee clobber
	make -C Zmumu clobber
	make -C ZmumuNLO clobber
	make -C Wenu clobber
	make -C Wmunu clobber
	make -C WmunuNLO clobber
	make -C ttbar clobber
	make -C Htautau clobber
	make -C Ztautau clobber
	make -C ScalNLO clobber
