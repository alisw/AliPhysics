# Top level Makefile for AliRoot System
#
# Author: Jan-Erik Revsbech (revsbech@fys.ku.dk)
#         Developed on idea of Boris Polichtchouk (Boris.Polichtchouk@cern.ch), 15/4/2001

# /* $Id$ */

ifdef ALIVERBOSE
MUTE:=
else
MUTE:=@
endif

CLEAN=$(findstring clean,$(MAKECMDGOALS))

#-------------------------------------------------------------------------------
# IRST coding rule check

IRST_INSTALLDIR=$(ALICE)/local/IRST
IRST_CONFIG_DIR=$(IRST_INSTALLDIR)/userConfig/ALICE
CLASSPATH=$(IRST_INSTALLDIR)
export CLASSPATH IRST_INSTALLDIR IRST_CONFIG_DIR
CODE_CHECK=java rules.ALICE.ALICERuleChecker
REV_ENG=$(IRST_INSTALLDIR)/scripts/revEng.sh

SMELL_DETECTOR_DIR=$(IRST_INSTALLDIR)/smell-detector


#-------------------------------------------------------------------------------
# Include machine dependent macros

-include build/Makefile.$(ALICE_TARGET)

#-------------------------------------------------------------------------------
# Check if called with debug

ifeq ($(ALIDEBUG),YES)
override ALICE_TARGET := $(ALICE_TARGET)DEBUG
FFLAGS   := -g $(filter-out -O%,$(FFLAGS))
CXXFLAGS := -g $(filter-out -O%,$(CXXFLAGS))
CFLAGS   := -g $(filter-out -O%,$(CFLAGS))
SOFLAGS  := -g $(filter-out -O%,$(SOFLAGS))
LDFLAGS  := -g $(filter-out -O%,$(LDFLAGS))
endif

#-------------------------------------------------------------------------------
# Check if called with profile

ifeq ($(ALIPROFILE),YES)
override ALICE_TARGET := $(ALICE_TARGET)PROF
FFLAGS   += -pg
CXXFLAGS += -pg
CFLAGS   += -pg
SOFLAGS  += -pg
LDFLAGS  += -pg
endif

#-------------------------------------------------------------------------------
# Check if DATE is installed

ifdef DATE_ROOT
DATEFLAGS  = -DALI_DATE -D${DATE_SYS} -DDATE_SYS=${DATE_SYS} -Dlong32=${DATE_LONG32} \
             -Dlong64='${DATE_LONG64}' -DdatePointer=${DATE_POINTER} -I${DATE_COMMON_DEFS}
CXXFLAGS  += $(DATEFLAGS)
CFLAGS    += $(DATEFLAGS)
CINTFLAGS += $(DATEFLAGS)
DEPINC    += $(DATEFLAGS)
else
DATEFLAGS  = -D`uname` -DDATE_SYS=`uname` -Dlong32='int' \
             -Dlong64='long long' -DdatePointer='long'
CXXFLAGS  += $(DATEFLAGS)
CFLAGS    += $(DATEFLAGS)
CINTFLAGS += $(DATEFLAGS)
DEPINC    += $(DATEFLAGS)
endif

#-------------------------------------------------------------------------------
# Location where to install libraries and binaries and common header files

LIBPATH      := lib/tgt_$(ALICE_TARGET)
BINPATH      := bin/tgt_$(ALICE_TARGET)
EXPORTDIR    := $(ALICE_ROOT)/include
BINLIBDIRS   := -L$(ALICE_ROOT)/$(LIBPATH)

#-------------------------------------------------------------------------------
# Modules to build

ifeq (ALIVERBOSE,2)
$(warning MAKECMDGOALS=$(MAKECMDGOALS))
endif

ALIROOTMODULES := STEER PHOS TRD TPC ZDC MUON PMD FMD TOF ITS \
      CRT HMPID START STRUCT EVGEN RALICE VZERO \
      THijing MEVSIM TMEVSIM THbtp HBTP EMCAL HBTAN \
      THerwig TEPEMGEN EPEMGEN FASTSIM TPHIC RAW MONITOR ANALYSIS \
      JETAN HLT LHC EVE

ifeq ($(findstring TFluka,$(MAKECMDGOALS)),TFluka)
ALIROOTMODULES += TFluka
endif

ifeq ($(findstring PWG0,$(MAKECMDGOALS)),PWG0)
ALIROOTMODULES += PWG0
endif

ifeq ($(findstring PWG2,$(MAKECMDGOALS)),PWG2)
ALIROOTMODULES += PWG2
endif

ifeq ($(findstring PWG3,$(MAKECMDGOALS)),PWG3)
ALIROOTMODULES += PWG3
endif

ifeq ($(findstring SHUTTLE,$(MAKECMDGOALS)),SHUTTLE)
ALIROOTMODULES += SHUTTLE
endif

ifeq ($(findstring Flugg,$(MAKECMDGOALS)),Flugg)
ALIROOTMODULES += Flugg
endif

CERNMODULES := LHAPDF PYTHIA6 HIJING MICROCERN HERWIG

MODULES := $(ALIROOTMODULES) $(CERNMODULES) ALIROOT

MODDIRS := $(MODULES)

#-------------------------------------------------------------------------------
# Default include dirs for C++, Fortran, Cint, and dependencies
# The module directory will be added by each module

GENINC     := -I$(ALICE_ROOT)/include -I$(shell root-config --incdir)
CXXFLAGS   += $(GENINC)
CXXFLAGSNO += $(GENINC)
CINTFLAGS  += $(GENINC)
DEPINC     += $(GENINC)

#-------------------------------------------------------------------------------
# Libraries to link binaries against
# Libraries will be linked against SHLIB
# ROOT libraries 

ROOTCLIBS     := $(shell root-config --glibs) -lThread -lMinuit -lHtml -lVMC -lEG -lGeom -lTreePlayer

ROOTPLIBS     := -lEGPythia6

ALILIBS	      := -L$(LIBDIR) -lMUON -lTPC -lPMD -lTRD -lFMD -lTOF \
                -lITS -lPHOS -lCRT -lHMPID -lVZERO -lZDC -lSTRUCT \
                -lSTART -lEVGEN -lSTEER

LIBS := $(ROOTCLIBS) $(ROOTPLIBS) $(SYSLIBS)

#-------------------------------------------------------------------------------
# default target

default:
	$(MUTE)$(MAKE) aliroot

FORCE:

#-------------------------------------------------------------------------------
# Each module will add to these macros

ALLLIBS      :=
ALLEXECS     :=
INCLUDEFILES :=
BINLIBS      :=
EXPORTFILES  :=

#-------------------------------------------------------------------------------
# Dependencies of module.mk files if not cleaning

ifeq ($(CLEAN),)
include build/module.dep
endif

#-------------------------------------------------------------------------------
# Check if module.mk is present for the library

%.mk: build/module.tpl build/header.tpl build/clean.tpl share/alibtool
ifndef ALIQUIET
	@echo "***** Creating $@ file *****";
endif
	@share/alibtool mkmodule  $(patsubst %/module.mk,%,$@) > $@;

#-------------------------------------------------------------------------------
# If making modules, not not include anything

ifeq ($(findstring modules,$(MAKECMDGOALS)),)

#-------------------------------------------------------------------------------
# Include the modules

-include $(patsubst %,%/module.mk,$(MODULES))

#-------------------------------------------------------------------------------
# If cleaning, do not include dependencies or module.mk files.

ifeq ($(CLEAN),)

#-------------------------------------------------------------------------------
# Include dependencies if not making them!

ifneq ($(MAKECMDGOALS),depend)
ifneq ($(MAKECMDGOALS),)

ifeq (ALIVERBOSE,2)
$(warning INCLUDEFILES=$(INCLUDEFILES))
endif
-include $(INCLUDEFILES)

endif
endif
endif
endif

#-------------------------------------------------------------------------------
# Include dummy dependency file *MUST* be last includefile

include build/dummy.d


#-------------------------------------------------------------------------------
# Targets

.PHONY:		alilibs aliroot makedistr clean distclean clean-all \
		htmldoc profile modules depend

modules: $(patsubst %,%/module.mk,$(MODULES)) 

ifeq ($(ALIPROFILE),YES)
alilibs: $(LIBPATH) modules $(ALLLIBS) $(ALLALIBS)
else
alilibs: $(LIBPATH) modules $(ALLLIBS)
endif

aliroot: alilibs $(BINPATH) $(ALLEXECS) 

alimdc-static: $(LIBPATH) $(RAWDatabaseALIB) $(MDCALIB) $(ESDALIB)
	 $(MUTE)rm -rf $(LIBPATH)/libAliMDC.a
	 $(MUTE)rm -rf junk
	 mkdir junk && cd junk && ar x ../$(RAWDatabaseALIB) && ar x ../$(MDCALIB) && ar x ../$(ESDALIB) && ar r ../$(LIBPATH)/libAliMDC.a *.o && cd .. && rm -rf junk

#-------------------------------------------------------------------------------
# Single Makefile "distribution": Makefile + modules + mkdepend scripts
makedistr: $(MODULES)
	 tar -cvf MakeDistr.tar $(patsubst %,%/*.pkg,$(MODULES)) \
		Makefile create build/*

all: aliroot

depend: $(INCLUDEFILES)

debug:
ifndef ALIQUIET
	@echo "***** Entering DEBUG mode. *****"
endif
	@(export ALIDEBUG=YES && $(MAKE) aliroot)

profile:
ifndef ALIQUIET
	@echo "***** Entering PROFILE mode. *****"
endif
	@(export ALIPROFILE=YES && $(MAKE) aliroot)

$(MODULES):
ifndef ALIQUIET
	@echo "***** Making $@ *****"
endif
	@mkdir -p $@

$(BINPATH):
ifndef ALIQUIET
	@echo "***** Making $@ *****"
endif
	@mkdir -p $@

$(LIBPATH):
ifndef ALIQUIET
	@echo "***** Making $@ *****"
endif
	@mkdir -p $@

build/dummy.d: $(EXPORTFILES)
	@(if [ ! -f $@ ] ; then \
	   touch $@; \
	fi)

clean:
	@echo "***** No target clean, use one of these *****"
	@echo "	clean-aliroot     : Clean up all aliroot libraries"
	@echo "	clean-MODULENAME  : Clean everything from module MODULENAME"
	@echo "	clean-all         : Cleans up everything, including cern libraires"
	@echo "	distclean         : Like clean-all + clean all tgt_*'s"
	@echo "	clean-modules     : Clean all module.mk files in all modules"
	@echo "	clean-libs        : Clean all libraries (not object files)"
	@echo "********************************************"

clean-all: clean-modules clean-libs clean-bins
ifndef ALIQUIET
	@echo "***** Cleaning up everything ****"
endif
	$(MUTE)rm -rf $(patsubst %,%/tgt_$(ALICE_TARGET),$(MODULES))
	$(MUTE)rm -rf $(EXPORTDIR)

distclean: clean-all
	$(MUTE)rm -rf */tgt_* bin lib

#-------------------------------------------------------------------------------
# This cleans only libraries that are not CERN-libraries

clean-aliroot:   $(patsubst %,%/module.mk,$(ALIROOTMODULES)) $(patsubst %,clean-%,$(ALIROOTMODULES))

CHECKMODULES := $(MODULES)
CHECKMODULES := $(filter-out HBTP,$(CHECKMODULES))
CHECKMODULES := $(filter-out MEVSIM,$(CHECKMODULES))
CHECKMODULES := $(filter-out EPEMGEN,$(CHECKMODULES))
CHECKMODULES := $(filter-out TPHIC,$(CHECKMODULES))
CHECKMODULES := $(filter-out LHAPDF,$(CHECKMODULES))
CHECKMODULES := $(filter-out MICROCERN,$(CHECKMODULES))

check-all:    $(patsubst %,%/module.mk,$(CHECKMODULES)) $(patsubst %,check-%,$(CHECKMODULES))

smell-all:    $(patsubst %,%/module.mk,$(CHECKMODULES)) $(patsubst %,smell-%,$(CHECKMODULES))

reveng-all:   $(patsubst %,%/module.mk,$(CHECKMODULES)) $(patsubst %,reveng-%,$(CHECKMODULES))

revdisp-all:  $(patsubst %,%/module.mk,$(CHECKMODULES)) $(patsubst %,revdisp-%,$(CHECKMODULES))

smell-all:    $(patsubst %,%/module.mk,$(CHECKMODULES)) $(patsubst %,smell-%,$(CHECKMODULES))

clean-dicts:
ifndef ALIQUIET
	@echo "***** Cleaning up G__ files *****"
endif
	$(MUTE)rm -rf */tgt_$(ALICE_TARGET)/G__*

clean-modules:
ifndef ALIQUIET
	@echo "***** Cleaning up module.mk files *****"
endif
	$(MUTE)rm -rf $(patsubst %,%/module.mk,$(MODULES))

clean-depend:
ifndef ALIQUIET
	@echo "***** Cleaning up dependencies *****"
endif
	$(MUTE)echo rm `find . -name "*.d"`

clean-objects:
ifndef ALIQUIET
	@echo "***** Cleaning up .o files *****"
endif
	$(MUTE)echo rm `find . -name "*.o"`

clean-libs:
ifndef ALIQUIET
	@echo "***** Cleaning up library files *****"
endif
	$(MUTE)rm -rf lib/tgt_$(ALICE_TARGET)/*

clean-bins:
ifndef ALIQUIET
	@echo "***** Cleaning up binary files *****"
endif
	$(MUTE)rm -rf bin/tgt_$(ALICE_TARGET)

clean-check-all:  $(patsubst %,%/module.mk,$(CHECKMODULES)) $(patsubst %,clean-check-%,$(CHECKMODULES))

clean-smell-all:  $(patsubst %,%/module.mk,$(CHECKMODULES)) $(patsubst %,clean-smell-%,$(CHECKMODULES))

clean-reveng-all: $(patsubst %,%/module.mk,$(CHECKMODULES)) $(patsubst %,clean-reveng-%,$(CHECKMODULES))

htmldoc:
	@rm -rf html/roothtml
	@rm -f  html/picts
	@rm -f /tmp/macros
	@cd html ;\
	aliroot -q -b "mkhtml.C(0,1)" ;\
	ls ../macros/*.C > /tmp/macros ;\
	for i in $(ALIROOTMODULES) ; do \
		ls ../$$i/*.C 2>/dev/null >> /tmp/macros ;\
	done ;\
	for i in `cat /tmp/macros` ; do \
		echo $$i ; \
		aliroot -b -q "mkhtml.C(\"$$i\")" > /dev/null ;\
	done ;\
	./makeExampleList ;
	@ln -s ../picts html/picts
	@ln -s ../../picts html/roothtml/picts
	@ln -s ../../../picts html/roothtml/src/picts
	@ln -s ../../../picts html/roothtml/examples/picts
