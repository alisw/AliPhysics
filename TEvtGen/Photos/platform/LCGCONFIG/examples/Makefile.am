exampledir = $(top_srcdir)/examples

noinst_LTLIBRARIES = libphotos_hepevt_example_interface.la
libphotos_hepevt_example_interface_la_SOURCES = photos_hepevt_example_interface.cxx

### HEPEVT example ###
example_PROGRAMS = photos_hepevt_example.exe

photos_hepevt_example_exe_SOURCES      = photos_hepevt_example.f
photos_hepevt_example_exe_FFLAGS       = -fno-automatic -fno-backslash -ffixed-line-length-132

INCLUDES = -I$(prefix)/include

AM_LDFLAGS = -R $(prefix)/lib

LDADD =                          \
$(FLIBS)                         \
$(prefix)/lib/libPhotosCxxInterface.so \
$(prefix)/lib/libPhotosFortran.so

photos_hepevt_example_exe_LDADD = -lphotos_hepevt_example_interface $(LDADD)

### All other examples require HepMC ###
if HAS_HEPMC

AM_LDFLAGS += -R $(HEPMC_DIR)/lib

INCLUDES += -I$(HEPMC_DIR)/include
LDADD    += -L$(HEPMC_DIR)/lib -lHepMC

photos_hepevt_example_exe_LDADD += -L$(HEPMC_DIR)/lib -lHepMC

example_PROGRAMS += photos_standalone_example.exe

photos_standalone_example_exe_SOURCES  = photos_standalone_example.cxx

### Photos gun example (require Pythia8) ###
if HAS_PYTHIA8

AM_LDFLAGS  += -R $(PYTHIA8_DIR)/lib/archive

INCLUDES += -I$(PYTHIA8_DIR)/include
LDADD    += -L$(PYTHIA8_DIR)/lib/archive -lpythia8 -llhapdfdummy -lhepmcinterface

example_PROGRAMS += single_photos_gun_example.exe  photosLCG_pythia_example.exe

single_photos_gun_example_exe_SOURCES  = single_photos_gun_example.cxx
photosLCG_pythia_example_exe_SOURCES   = photosLCG_pythia_example.cxx

### Pythia example (requires Pythia8 and MC-Tester) ###
if HAS_MCTESTER

AM_LDFLAGS  += -R $(MCTESTER_DIR)/lib

INCLUDES += -I$(MCTESTER_DIR)/include -I$(ROOTINC)
LDADD    += -L$(MCTESTER_DIR)/lib -lHEPEvent -lHepMCEvent -lMCTester $(ROOTLIBS)

example_PROGRAMS += photos_pythia_example.exe

photos_pythia_example_exe_SOURCES  = photos_pythia_example.cxx

if HAS_TAUOLA

AM_LDFLAGS  += -R $(TAUOLA_DIR)/lib

INCLUDES += -I$(TAUOLA_DIR)/include
LDADD    += -L$(TAUOLA_DIR)/lib -lTauolaFortran -lTauolaCxxInterface

example_PROGRAMS += tauola_photos_pythia_example.exe

tauola_photos_pythia_example_exe_SOURCES  = tauola_photos_pythia_example.cxx

endif # TAUOLA

endif # MC-TESTER

endif # PYTHIA8

endif # HEPMC
