# Makefile is a part of the PYTHIA event generator.
# Copyright (C) 2018 Torbjorn Sjostrand.
# PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
# Please respect the MCnet Guidelines, see GUIDELINES for details.
# Author: Philip Ilten, October 2014 - November 2017.
#
# This is is the Makefile used to build PYTHIA on POSIX systems.
# Example usage is:
#     make -j2
# For help using the make command please consult the local system documentation,
# i.e. "man make" or "make --help".

################################################################################
# VARIABLES: Definition of the relevant variables from the configuration script
# and the distribution structure.
################################################################################

# Set the shell.
SHELL=/usr/bin/env bash

# Include the configuration and set the local directory structure.
ifeq (,$(findstring clean, $(MAKECMDGOALS)))
  -include Makefile.inc
endif
LOCAL_BIN=bin
LOCAL_DOCS=AUTHORS COPYING GUIDELINES README ../../examples/Makefile.inc
LOCAL_EXAMPLE=examples
LOCAL_INCLUDE=include
LOCAL_LIB=lib
LOCAL_SHARE=share/Pythia8
LOCAL_SRC=src
LOCAL_TMP=tmp
LOCAL_MKDIRS:=$(shell mkdir -p $(LOCAL_TMP) $(LOCAL_LIB))
CXX_COMMON:=-I$(LOCAL_INCLUDE) $(CXX_COMMON)

# PYTHIA.
OBJECTS=$(patsubst $(LOCAL_SRC)/%.cc,$(LOCAL_TMP)/%.o,\
	$(sort $(wildcard $(LOCAL_SRC)/*.cc)))
TARGETS=$(LOCAL_LIB)/libpythia8.a
ifeq ($(ENABLE_SHARED),true)
  TARGETS+=$(LOCAL_LIB)/libpythia8$(LIB_SUFFIX)
endif

# LHAPDF.
CXX_LHAPDF:=$(CXX_COMMON)
ifeq ($(LHAPDF5_USE),true)
  CXX_LHAPDF:=-I$(LHAPDF5_INCLUDE) $(CXX_LHAPDF)
  TARGETS+=$(LOCAL_LIB)/libpythia8lhapdf5.so
endif
ifeq ($(LHAPDF6_USE),true)
  CXX_LHAPDF:=-I$(LHAPDF6_INCLUDE) $(CXX_LHAPDF)
  TARGETS+=$(LOCAL_LIB)/libpythia8lhapdf6.so
endif
ifeq ($(BOOST_USE),true)
  CXX_LHAPDF:=-I$(BOOST_INCLUDE) $(CXX_LHAPDF)
endif

# POWHEG (needs directory that contains just POWHEG binaries and scripts).
ifeq ($(POWHEG_USE),true)
  ifneq ($(POWHEG_DIR),./)
    TARGETS+=$(patsubst $(POWHEG_BIN)/%,$(LOCAL_LIB)/libpythia8powheg%.so,\
             $(wildcard $(POWHEG_BIN)/*))
  endif
endif

# Python.
PYTHON_COMMON=-I$(PYTHON_INCLUDE) $(CXX_COMMON) -Wl,-rpath,$(PREFIX_LIB)
ifeq ($(PYTHON_USE),true)
  TARGETS+=$(LOCAL_LIB)/_pythia8.so
endif

# GZIP.
OBJ_COMMON=-MD $(CXX_COMMON)
LIB_COMMON=-Wl,-rpath,$(PREFIX_LIB) -ldl
ifeq ($(GZIP_USE),true)
  PYTHON_COMMON+= -DGZIPSUPPORT -I$(GZIP_INCLUDE)
  PYTHON_COMMON+= -L$(GZIP_LIB) -Wl,-rpath,$(GZIP_LIB) -lz
  OBJ_COMMON+= -DGZIPSUPPORT -I$(GZIP_INCLUDE)
  LIB_COMMON+= -L$(GZIP_LIB) -Wl,-rpath,$(GZIP_LIB) -lz
endif

################################################################################
# RULES: Definition of the rules used to build PYTHIA.
################################################################################

# Rules without physical targets (secondary expansion for documentation).
.SECONDEXPANSION:
.PHONY: all install clean distclean

# All targets.
all: $(TARGETS) $(addprefix $(LOCAL_SHARE)/, $(LOCAL_DOCS))

# The documentation.
$(addprefix $(LOCAL_SHARE)/, $(LOCAL_DOCS)): $$(notdir $$@)
	cp $^ $@

# The Makefile configuration.
Makefile.inc:
	./configure

# Auto-generated (with -MD flag) dependencies.
-include $(LOCAL_TMP)/*.d

# PYTHIA.
$(LOCAL_TMP)/Pythia.o: $(LOCAL_SRC)/Pythia.cc Makefile.inc
	$(CXX) $< -o $@ -c $(OBJ_COMMON) -DXMLDIR=\"$(PREFIX_SHARE)/xmldoc\"
$(LOCAL_TMP)/%.o: $(LOCAL_SRC)/%.cc
	$(CXX) $< -o $@ -c $(OBJ_COMMON)
$(LOCAL_LIB)/libpythia8.a: $(OBJECTS)
	rm -f $(LOCAL_LIB)/libpythia8$(LIB_SUFFIX)
	ar cru $@ $^
$(LOCAL_LIB)/libpythia8$(LIB_SUFFIX): $(OBJECTS)
	$(CXX) $^ -o $@ $(CXX_COMMON) $(CXX_SHARED) $(CXX_SONAME)$(notdir $@)\
	  $(LIB_COMMON)

# LHAPDF (turn off all warnings for readability).
$(LOCAL_TMP)/LHAPDF%Plugin.o: $(LOCAL_INCLUDE)/Pythia8Plugins/LHAPDF%.h
	$(CXX) -x c++ $< -o $@ -c -MD -w $(CXX_LHAPDF)
$(LOCAL_LIB)/libpythia8lhapdf5.so: $(LOCAL_TMP)/LHAPDF5Plugin.o\
	$(LOCAL_LIB)/libpythia8.a
	$(CXX) $^ -o $@ $(CXX_COMMON) $(CXX_SHARED) $(CXX_SONAME)$(notdir $@)\
	 -L$(LHAPDF5_LIB) -Wl,-rpath,$(LHAPDF5_LIB) -lLHAPDF -lgfortran
$(LOCAL_LIB)/libpythia8lhapdf6.so: $(LOCAL_TMP)/LHAPDF6Plugin.o\
	$(LOCAL_LIB)/libpythia8.a
	$(CXX) $^ -o $@ $(CXX_COMMON) $(CXX_SHARED) $(CXX_SONAME)$(notdir $@)\
	 -L$(LHAPDF6_LIB) -Wl,-rpath,$(LHAPDF6_LIB) -lLHAPDF

# POWHEG (exclude any executable ending with sh).
$(LOCAL_TMP)/POWHEGPlugin.o: $(LOCAL_INCLUDE)/Pythia8Plugins/LHAPowheg.h
	$(CXX) -x c++ $< -o $@ -c -MD -w $(CXX_COMMON)
$(LOCAL_LIB)/libpythia8powheg%sh.so: $(POWHEG_BIN)/%sh;
$(LOCAL_LIB)/libpythia8powheg%.so: $(POWHEG_BIN)/% $(LOCAL_TMP)/POWHEGPlugin.o\
	$(LOCAL_LIB)/libpythia8.a
	ln -s $< $(notdir $<); $(CXX) $(notdir $<) $(LOCAL_TMP)/POWHEGPlugin.o\
	 $(LOCAL_LIB)/libpythia8.a -o $@ $(CXX_COMMON) $(CXX_SHARED)\
	 $(CXX_SONAME)$(notdir $@) -Wl,-rpath,$(POWHEG_BIN); rm $(notdir $<)

# Python (turn off all warnings for readability).
$(LOCAL_LIB)/pythia8.py: $(LOCAL_INCLUDE)/Pythia8Plugins/PythonWrapper.h
	SPLIT=`grep -n "PYTHON SOURCE" $< | cut -d : -f 1`;\
	 SPLIT=$$[$$SPLIT+1]; tail -n +$$SPLIT $< | cut -d "/" -f 3- > $@
	$(PYTHON_BIN)python -m compileall $(LOCAL_LIB)
$(LOCAL_LIB)/_pythia8.so: $(LOCAL_INCLUDE)/Pythia8Plugins/PythonWrapper.h\
	$(LOCAL_LIB)/pythia8.py $(wildcard $(LOCAL_INCLUDE)/*/*.h) |\
	$(LOCAL_LIB)/libpythia8$(LIB_SUFFIX)
	$(CXX) -x c++ $< -o $@ -w $(PYTHON_COMMON) $(CXX_SHARED)\
	 -Wl,-undefined,dynamic_lookup -Wno-long-long\
	 $(CXX_SONAME)$(notdir $@) -L$(LOCAL_LIB) -lpythia8

# Install (rsync is used for finer control).
install: all
	mkdir -p $(PREFIX_BIN) $(PREFIX_INCLUDE) $(PREFIX_LIB) $(PREFIX_SHARE)
	rm -f $(PREFIX_LIB)/libpythia8$(LIB_SUFFIX)
	rsync -a $(LOCAL_BIN)/* $(PREFIX_BIN) --exclude .svn
	rsync -a $(LOCAL_INCLUDE)/* $(PREFIX_INCLUDE) --exclude .svn
	rsync -a $(LOCAL_LIB)/* $(PREFIX_LIB) --exclude .svn
	rsync -a $(LOCAL_SHARE)/* $(PREFIX_SHARE) --exclude .svn
	rsync -a $(LOCAL_EXAMPLE) $(PREFIX_SHARE) --exclude .svn

# Clean.
clean:
	rm -rf $(LOCAL_TMP) $(LOCAL_LIB)
	rm -f $(LOCAL_EXAMPLE)/*Dct.*
	rm -f $(LOCAL_EXAMPLE)/*[0-9]
	rm -f $(LOCAL_EXAMPLE)/weakbosons.lhe
	rm -f $(LOCAL_EXAMPLE)/Pythia8.promc
	rm -f $(LOCAL_EXAMPLE)/hist.root

# Clean all temporary and generated files.
distclean: clean
	find . -type f -name Makefile.inc -print0 | xargs -0 rm -f
	find . -type f -name "*~" -print0 | xargs -0 rm -f
	find . -type f -name "#*" -print0 | xargs -0 rm -f
	rm -f $(LOCAL_SHARE)/AUTHORS
	rm -f $(LOCAL_SHARE)/COPYING
	rm -f $(LOCAL_SHARE)/GUIDELINES
	rm -f $(LOCAL_SHARE)/README
