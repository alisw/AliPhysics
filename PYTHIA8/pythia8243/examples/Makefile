# Makefile is a part of the PYTHIA event generator.
# Copyright (C) 2018 Torbjorn Sjostrand.
# PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
# Please respect the MCnet Guidelines, see GUIDELINES for details.
# Author: Philip Ilten, September 2014.
#
# This is is the Makefile used to build PYTHIA examples on POSIX systems.
# Example usage is:
#     make main01
# For help using the make command please consult the local system documentation,
# i.e. "man make" or "make --help".

################################################################################
# VARIABLES: Definition of the relevant variables from the configuration script.
################################################################################

# Set the shell.
SHELL=/usr/bin/env bash

# Include the configuration.
-include Makefile.inc

# Handle GZIP support.
GZIP_INC=
GZIP_FLAGS=
ifeq ($(GZIP_USE),true)
  GZIP_INC+= -DGZIPSUPPORT -I$(GZIP_INCLUDE)
  GZIP_FLAGS+= -L$(GZIP_LIB) -Wl,-rpath,$(GZIP_LIB) -lz
endif

# Check distribution (use local version first, then installed version).
ifneq ("$(wildcard ../lib/libpythia8.*)","")
  PREFIX_LIB=../lib
  PREFIX_INCLUDE=../include
endif
CXX_COMMON:=-I$(PREFIX_INCLUDE) $(CXX_COMMON)
CXX_COMMON+= -L$(PREFIX_LIB) -Wl,-rpath,$(PREFIX_LIB) -lpythia8 -ldl 

################################################################################
# RULES: Definition of the rules used to build the PYTHIA examples.
################################################################################

# Rules without physical targets (secondary expansion for specific rules).
.SECONDEXPANSION:
.PHONY: all clean

# All targets (no default behavior).
all:
	@echo "Usage: make mainXX"

# The Makefile configuration.
Makefile.inc:
	$(error Error: PYTHIA must be configured, please run "./configure"\
                in the top PYTHIA directory)

# PYTHIA libraries.
$(PREFIX_LIB)/libpythia8.a :
	$(error Error: PYTHIA must be built, please run "make"\
                in the top PYTHIA directory)

# Examples without external dependencies.
main% : main%.cc $(PREFIX_LIB)/libpythia8.a
	$(CXX) $< -o $@ $(CXX_COMMON) $(GZIP_INC) $(GZIP_FLAGS)

# MixMax.
main23: $$@.cc $(PREFIX_LIB)/libpythia8.a
	$(CXX) $< -o $@ -std=c++11 -w $(CXX_COMMON) $(GZIP_INC) $(GZIP_FLAGS)

# GZIP (required).
main34: $$@.cc $(PREFIX_LIB)/libpythia8.a
ifeq ($(GZIP_USE),true)
	$(CXX) $< -o $@ $(CXX_COMMON) $(GZIP_INC) $(GZIP_FLAGS)
else
	@echo "Error: $@ requires GZIP"
endif

# HEPMC2.
main41 main42 main43 main85 main86 main87 main88 main89: $$@.cc\
	$(PREFIX_LIB)/libpythia8.a
ifeq ($(HEPMC2_USE),true)
	$(CXX) $< -o $@ -I$(HEPMC2_INCLUDE) $(CXX_COMMON)\
	 -L$(HEPMC2_LIB) -Wl,-rpath,$(HEPMC2_LIB) -lHepMC\
	 $(GZIP_INC) $(GZIP_FLAGS)
else
	@echo "Error: $@ requires HEPMC2"
endif

# PROMC.
main46: $$@.cc $(PREFIX_LIB)/libpythia8.a
ifeq ($(PROMC_USE),true)
	$(CXX) $< -o $@ -I$(PROMC_INCLUDE)/src -I$(PROMC_INCLUDE)/include\
	 $(CXX_COMMON) -DPROMC=\"$(PROMC_INCLUDE)\" -Wno-long-long\
	 -L$(PROMC_LIB) -Wl,-rpath,$(PROMC_LIB) -lpromc -lprotoc -lprotobuf\
	 -lprotobuf-lite -lcbook $(GZIP_INC) $(GZIP_FLAGS)
else
	@echo "Error: $@ requires PROMC"
endif

# EVTGEN (and HEPMC2).
main48: $$@.cc $(PREFIX_LIB)/libpythia8.a
ifeq ($(EVTGEN_USE)$(HEPMC2_USE)$(ENABLE_SHARED),truetruetrue)
	$(CXX) $< -o $@ -I$(EVTGEN_INCLUDE) $(CXX_COMMON)\
	 -DEVTGEN_PYTHIA -DEVTGEN_EXTERNAL -Wl,-rpath,$(HEPMC2_LIB)\
	 -L$(EVTGEN_LIB) -Wl,-rpath,$(EVTGEN_LIB) -lEvtGenExternal -lEvtGen\
	  $(GZIP_INC) $(GZIP_FLAGS)
else
	@echo "Error: $@ requires EVTGEN and HEPMC2"
endif

# FASTJET3.
main71 main72 main75: $$@.cc $(PREFIX_LIB)/libpythia8.a
ifeq ($(FASTJET3_USE),true)
	$(CXX) $< -o $@ -I$(FASTJET3_INCLUDE) $(CXX_COMMON)\
	 -L$(FASTJET3_LIB) -Wl,-rpath,$(FASTJET3_LIB) -lfastjet\
	  $(GZIP_INC) $(GZIP_FLAGS)
else
	@echo "Error: $@ requires FASTJET3"
endif

# FASTJET3 with modified Mass-Drop Tagger.
main74: $$@.cc $(PREFIX_LIB)/libpythia8.a
ifeq ($(FASTJET3_USE),true)
	$(CXX) $< -o $@ -I$(FASTJET3_INCLUDE) $(CXX_COMMON)\
	 -L$(FASTJET3_LIB) -Wl,-rpath,$(FASTJET3_LIB) -lfastjet -lRecursiveTools\
	  $(GZIP_INC) $(GZIP_FLAGS)
else
	@echo "Error: $@ requires FASTJET3"
endif

# FASTJET3 and HEPMC2.
main81 main82 main83 main84: $$@.cc $(PREFIX_LIB)/libpythia8.a
ifeq ($(FASTJET3_USE)$(HEPMC2_USE),truetrue)
	$(CXX) $< -o $@ -I$(FASTJET3_INCLUDE) -I$(HEPMC2_INCLUDE) $(CXX_COMMON)\
	 -L$(HEPMC2_LIB) -Wl,-rpath,$(HEPMC2_LIB) -lHepMC\
	 -L$(FASTJET3_LIB) -Wl,-rpath,$(FASTJET3_LIB) -lfastjet\
	  $(GZIP_INC) $(GZIP_FLAGS)
else
	@echo "Error: $@ requires FASTJET3 and HEPMC2"
endif

# ROOT (turn off all warnings for readability).
main91: $$@.cc $(PREFIX_LIB)/libpythia8.a
ifeq ($(ROOT_USE),true)
	$(CXX) $< -o $@ -w -I$(ROOT_INCLUDE) $(CXX_COMMON)\
	 `$(ROOTBIN)root-config --cflags`\
	 -Wl,-rpath,$(ROOT_LIB) `$(ROOT_BIN)root-config --glibs`
else
	@echo "Error: $@ requires ROOT"
endif
main92: $$@.cc $(PREFIX_LIB)/libpythia8.a main92.so
ifeq ($(ROOT_USE),true)
	$(CXX) $< main92.so -o $@ -w -I$(ROOT_INCLUDE) $(CXX_COMMON)\
	 `$(ROOTBIN)root-config --cflags` -Wl,-rpath,./\
	 -Wl,-rpath,$(ROOT_LIB) `$(ROOT_BIN)root-config --glibs`
else
	@echo "Error: $@ requires ROOT"
endif
main92.so: main92Dct.cc $(PREFIX_LIB)/libpythia8.a
	$(CXX) $< -o $@ -w -I$(ROOT_INCLUDE) $(CXX_SHARED) $(CXX_COMMON)\
	 `$(ROOTBIN)root-config --cflags`
main92Dct.cc: main92.h main92LinkDef.h
	export LD_LIBRARY_PATH=$$LD_LIBRARY_PATH:$(ROOT_LIB);\
	 $(ROOT_BIN)rootcint -f $@ -c -I$(PREFIX_INCLUDE) $^

# main93 with several dependencies. It is assumed that
# Rivet and YODA are installed to system path or appended to it.
main93: $$@.cc $(PREFIX_LIB)/libpythia8.a main93.so
ifeq ($(ROOT_USE),true)
	$(CXX) $< main93.so -o $@ -DUSE_ROOT -DUSE_YODA -std=c++11 -w -I$(ROOT_INCLUDE) $(CXX_COMMON)\
	 `$(ROOTBIN)root-config --cflags` -Wl,-rpath,./\
	 -Wl,-rpath,$(ROOT_LIB) `$(ROOT_BIN)root-config --glibs` -lYODA\
	 -lHepMC -lRivet  $(GZIP_INC) $(GZIP_FLAGS)
else
	$(CXX) $< -o $@ -DUSE_YODA -std=c++11 -w $(CXX_COMMON) -lYODA\
	 -lHepMC -lRivet  $(GZIP_INC) $(GZIP_FLAGS)
endif
main93.so: main93Dct.cc $(PREFIX_LIB)/libpythia8.a
ifeq ($(ROOT_USE),true)
	$(CXX) $< -o $@ -w -I$(ROOT_INCLUDE) $(CXX_SHARED) $(CXX_COMMON)\
	 `$(ROOTBIN)root-config --cflags`
else
	@echo "Skipping Root shared library"
endif
main93Dct.cc: main93.h main93LinkDef.h
ifeq ($(ROOT_USE),true)
	export LD_LIBRARY_PATH=$$LD_LIBRARY_PATH:$(ROOT_LIB);\
	 $(ROOT_BIN)rootcint -f $@ -c -I$(PREFIX_INCLUDE) $^
else
	@echo "Skipping Root dictionary"
endif

# User-written examples for tutorials, without external dependencies.
mymain% : mymain%.cc $(PREFIX_LIB)/libpythia8.a
	$(CXX) $< -o $@ $(CXX_COMMON) $(GZIP_INC) $(GZIP_FLAGS)

# Internally used tests, without external dependencies.
test% : test%.cc $(PREFIX_LIB)/libpythia8.a
	$(CXX) $< -o $@ $(CXX_COMMON) $(GZIP_INC) $(GZIP_FLAGS)

# Clean.
clean:
	@rm -f main[0-9][0-9]; rm -f out[0-9][0-9];\
	rm -f main[0-9][0-9][0-9]; rm -f out[0-9][0-9][0-9];\
	rm -f mymain[0-9][0-9]; rm -f myout[0-9][0-9];\
	rm -f test[0-9][0-9][0-9]; rm -f *.dat;\
	rm -f weakbosons.lhe; rm -f Pythia8.promc; rm -f hist.root;\
	rm -f *~; rm -f \#*; rm -f core*; rm -f *Dct.*; rm -f *.so;
