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
# The compilers
CXX           := $(shell root-config --cxx)
F77	      := $(shell root-config --f77)
CC	      := $(shell root-config --cc)

CXXWARN       =

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

ifeq ($(shell which date-config 2>/dev/null),)
DATEFLAGS := -D$(shell uname) -DDATE_SYS=$(shell uname) -Dlong32='int' \
             -Dlong64='long long' -DdatePointer='long'
DMONLIBS  :=
else 
DATEFLAGS := -DALI_DATE $(shell date-config --cflags | tr \" \')
DMONLIBS  := $(shell date-config --monitorlibs)
endif
CXXFLAGS  += $(DATEFLAGS)
CFLAGS    += $(DATEFLAGS)
CINTFLAGS += $(DATEFLAGS)
DEPINC    += $(DATEFLAGS)

#-------------------------------------------------------------------------------
# Add warning flags if any

ifneq ($(ALICXXWARN),no)
CXXFLAGS  += $(CXXWARN)
endif

#-------------------------------------------------------------------------------
# ROOT Stuff
ROOTCONFIG    := root-config
ROOTLIBDIR    := $(shell $(ROOTCONFIG) --libdir)
ROOTPLUGDIR   := $(ROOTLIBDIR)/$(dir $(shell $(ROOTCONFIG) --version))
ROOTINCDIR    := $(shell $(ROOTCONFIG) --incdir)
ROOTCLIBS     := $(shell $(ROOTCONFIG) --glibs) \
			-lThread 		\
			-lMinuit 		\
			-lVMC 			\
			-lEG 			\
			-lGeom 			\
			-lTreePlayer 		\
			-lXMLIO 		\
			-lXMLParser 		\
			-lProof 		\
			-lProofPlayer 		\
			-lMLP 			\
			-lSpectrum 		\
			-L$(ROOTPLUGDIR)
CHECKALIEN    := $(shell root-config --has-alien)
CHECKXML      := $(shell root-config --has-xml)

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
      ACORDE HMPID T0 BCM STRUCT EVGEN RALICE VZERO \
      THijing THbtp EMCAL \
      THerwig TEPEMGEN FASTSIM TPHIC RAW MONITOR ANALYSIS \
      JETAN HLT LHC ESDCheck STAT TTherminator CORRFW DPMJET TDPMjet \
      PWG0 PWG1 PWG2 PWG3 PWG4 TRIGGER OADB

# Additional generators
ALIROOTMODULES += TUHKMgen
ALIROOTMODULES += EPOS
ALIROOTMODULES += PYTHIA8
ALIROOTMODULES += TAmpt

ifneq ($(shell $(ROOTCONFIG) --has-opengl), no)
ALIROOTMODULES += EVE
endif 

ifeq ($(findstring TFluka,$(MAKECMDGOALS)),TFluka)
ALIROOTMODULES += TFluka
endif

ifeq ($(findstring THydjet,$(MAKECMDGOALS)),THydjet)
ALIROOTMODULES += THydjet
endif

ifeq ($(findstring SHUTTLE,$(MAKECMDGOALS)),SHUTTLE)
ALIROOTMODULES += SHUTTLE
endif

CERNMODULES := LHAPDF HIJING MICROCERN HERWIG	 
ifneq ($(wildcard $(ROOTINCDIR)/TPythia6.h),)
CERNMODULES += PYTHIA6
endif 

MODULES := $(ALIROOTMODULES) $(CERNMODULES) ALIROOT

MODDIRS := $(MODULES)

#-------------------------------------------------------------------------------
# Default include dirs for C++, Fortran, Cint, and dependencies
# The module directory will be added by each module

GENINC     := -I$(ALICE_ROOT)/include -isystem$(shell root-config --incdir)
RCFLAGS    := $(shell root-config --auxcflags) 
RLFLAGS    := $(shell root-config --ldflags)
CXXFLAGS   += $(GENINC) $(RCFLAGS)
CXXFLAGSNO += $(GENINC) $(RCFLAGS) -Wno-write-strings
CFLAGS     += $(GENINC) $(RCFLAGS)
CINTFLAGS  += $(GENINC) $(RCFLAGS)
FFLAGS	   += $(RCFLAGS)
ifeq (macosxicc,$(ALICE_TARGET))
FFLAGS     := $(patsubst -pthread, -reentrancy threaded, $(FFLAGS))
endif
LDFLAGS    += $(RLFLAGS)
SOFLAGS    += $(RLFLAGS)
DEPINC     += $(GENINC)

#-------------------------------------------------------------------------------
# Libraries to link binaries against
# Libraries will be linked against SHLIB
# ROOT libraries 

ALILIBS	      := -L$(LIBDIR) -lMUON -lTPC -lPMD -lTRD -lFMD -lTOF \
                -lITS -lPHOS -lACORDE -lHMPID -lVZERO -lZDC -lSTRUCT \
                -lT0 -lEVGEN -lSTEER -lTRIGGER

LIBS := $(ROOTCLIBS) $(ROOTPLIBS) $(SYSLIBS)

ARVERSIONFILE := $(EXPORTDIR)/ARVersion.h
SVNREV        := $(strip $(shell LANG=C LANGUAGE=C svn info | grep "Last Changed Rev:" | \
			         cut -d: -f2 ))
SVNBRANCH     := $(subst //alisoft.cern.ch/AliRoot/,,$(shell svn info | grep "URL:" | cut -d: -f3 ))

#-------------------------------------------------------------------------------
# default target

default: $(ARVERSIONFILE)
	$(MUTE)$(MAKE) aliroot

#-------------------------------------------------------------------------------
# Write header file with aliroot svn version and url

$(ARVERSIONFILE): $(ALICE_ROOT)/.svn/entries $(EXPORTDIR)
	$(MUTE)rm -f $(ARVERSIONFILE)
	@echo "***** Making $(ARVERSIONFILE) *****"
	@echo "#ifndef ALIROOT_ARVersion" >> $@
	@echo "#define ALIROOT_ARVersion" >> $@
	@echo "#define ALIROOT_SVN_REVISION $(SVNREV)" >> $@
	@echo "#define ALIROOT_SVN_BRANCH \"$(SVNBRANCH)\"" >> $@
	@echo "#endif" >> $@ 
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

ALIRECO.par: macros/loadlibsrec.C STEER/PROOF-INF.ALIRECO/SETUP.C
	$(MUTE)echo "***** Creating package archive" $@ "*****"
	$(MUTE)rm -rf ALIRECO
	$(MUTE)mkdir -p ALIRECO/PROOF-INF
	$(MUTE)cat $^ > ALIRECO/PROOF-INF/SETUP.C
	(tar cfzh $@ ALIRECO 2> /dev/null && echo "Package archive" $@ "created in" $(PWD)/$@) || (tar cfzh /tmp/$@ ALIRECO 2> /dev/null && echo "Package archive" $@ "created in /tmp/"$@)
	$(MUTE)rm -rf ALIRECO

ROOTALIBDIR:=$(shell root-config --libdir)

ALIMDCSPECFILE=$(RAWDIRO)/alimdc.spec
ALIMDCVERSION=$(subst -,.,$(notdir $(subst /RAW/mdc.h,,$(shell svn info RAW/mdc.h | grep "URL:" | cut -d: -f3 ))))
ALIMDCRELEASE=$(firstword $(shell svn info RAW/mdc.h | grep "Revision:" | cut -d: -f2 ))

alimdc-rpm: alimdc-static alimdc-specfile
	$(MUTE)rm -rf alimdc-root
	$(MUTE)mkdir -p alimdc-root/opt/alimdc/lib
	$(MUTE)mkdir -p alimdc-root/opt/alimdc/include
	$(MUTE)cp RAW/mdc.h alimdc-root/opt/alimdc/include
	$(MUTE)cp $(LIBPATH)/libAliMDC.a \
	$(ROOTALIBDIR)/libRoot.a \
	$(ROOTALIBDIR)/libfreetype.a $(ROOTALIBDIR)/libpcre.a \
	alimdc-root/opt/alimdc/lib
	$(MUTE)rm -rf RPMS
	$(MUTE)case `uname -m` in \
	    i?86*)	ALIMDCARCHDIR=i386;;\
	    ia64*)	ALIMDCARCHDIR=ia64;;\
	    x86_64*)	ALIMDCARCHDIR=x86_64;;\
	    *)		echo "Unknown architecture: `uname -m`"; exit 1;;\
	esac; \
	mkdir -p RPMS/$$ALIMDCARCHDIR; \
	rpmbuild --verbose --define "_topdir $(ALICE_ROOT)" --define "_tmppath $(ALICE_ROOT)" -bb $(ALIMDCSPECFILE); \
	cp -p RPMS/$$ALIMDCARCHDIR/alimdc-*.rpm .;
	$(MUTE)rm -rf alimdc-root
	$(MUTE)rm -rf RPMS
	@echo "***** alimdc RPM created and put $(ALICE_ROOT) folder *****"

alimdc-specfile: $(RAWDIRO)
	$(MUTE)rm -rf $(ALIMDCSPECFILE)
	@echo "***** Making alimdc RPM spec-file $(ALIMDCSPECFILE) *****"
	@echo "# RPM specfile for alimdc static libs" >> $(ALIMDCSPECFILE)
	@echo "# Package contains both ROOT and AliRoot" >> $(ALIMDCSPECFILE)
	@echo "# static libs needed by mStreamRecorder" >> $(ALIMDCSPECFILE)
	@echo "# in order to ROOT-ify the incoming raw" >> $(ALIMDCSPECFILE)
	@echo "# data" >> $(ALIMDCSPECFILE)
	@echo "# Example how-to build alimdc RPM:" >> $(ALIMDCSPECFILE)
	@echo "# cd $ALICE_ROOT" >> $(ALIMDCSPECFILE)
	@echo "# make alimdc-rpm" >> $(ALIMDCSPECFILE)
	@echo "" >> $(ALIMDCSPECFILE)
	@echo "Summary: AliMDC static libraries" >> $(ALIMDCSPECFILE)
	@echo "Name: alimdc" >> $(ALIMDCSPECFILE)
	@echo "Version:  $(ALIMDCVERSION)" >> $(ALIMDCSPECFILE)
	@echo "Release: $(ALIMDCRELEASE)" >> $(ALIMDCSPECFILE)
	@echo "# Copyright: CERN Alice Off-line" >> $(ALIMDCSPECFILE)
	@echo "License: CERN Alice Off-line" >> $(ALIMDCSPECFILE)
	@echo "Vendor: ALICE Core Off-line Group" >> $(ALIMDCSPECFILE)
	@echo "URL: http://aliceinfo.cern.ch" >> $(ALIMDCSPECFILE)
	@echo "Group: Applications/Alice" >> $(ALIMDCSPECFILE)
	@echo "Prefix: /opt/%{name}" >> $(ALIMDCSPECFILE)
	@echo "BuildRoot: %{_tmppath}/%{name}-root" >> $(ALIMDCSPECFILE)
	@echo "" >> $(ALIMDCSPECFILE)
	@echo "# automatic dependencies" >> $(ALIMDCSPECFILE)
	@echo "AutoReqProv: yes" >> $(ALIMDCSPECFILE)
	@echo "" >> $(ALIMDCSPECFILE)
	@echo "# list here required RPM packages for runtime" >> $(ALIMDCSPECFILE)
	@echo "Requires: glibc" >> $(ALIMDCSPECFILE)
	@echo "" >> $(ALIMDCSPECFILE)
	@echo "Provides: alimdc" >> $(ALIMDCSPECFILE)
	@echo "" >> $(ALIMDCSPECFILE)
	@echo "# description of the package" >> $(ALIMDCSPECFILE)
	@echo "%description" >> $(ALIMDCSPECFILE)
	@echo "Package contains both ROOT and AliRoot" >> $(ALIMDCSPECFILE)
	@echo "static libs needed by mStreamRecorder" >> $(ALIMDCSPECFILE)
	@echo "in order to ROOT-ify the incoming raw" >> $(ALIMDCSPECFILE)
	@echo "data. The package version correspond to" >> $(ALIMDCSPECFILE)
	@echo "the AliRoot one." >> $(ALIMDCSPECFILE)
	@echo "" >> $(ALIMDCSPECFILE)
	@echo "# list of files to be installed" >> $(ALIMDCSPECFILE)
	@echo "%files" >> $(ALIMDCSPECFILE)
	@echo "%defattr (-,root,root)" >> $(ALIMDCSPECFILE)
	@echo "%{prefix}/lib/libAliMDC.a" >> $(ALIMDCSPECFILE)
	@echo "%{prefix}/lib/libRoot.a" >> $(ALIMDCSPECFILE)
	@echo "%{prefix}/lib/libpcre.a" >> $(ALIMDCSPECFILE)
	@echo "%{prefix}/lib/libfreetype.a" >> $(ALIMDCSPECFILE)
	@echo "%{prefix}/include/mdc.h" >> $(ALIMDCSPECFILE)

alimdc-static: $(LIBPATH) $(BINPATH) $(RAWDatabaseALIB) $(MDCALIB) $(ESDALIB) $(STEERBaseALIB) $(alimdcCXXO)
	 $(MUTE)rm -rf $(LIBPATH)/libAliMDC.a
	 $(MUTE)rm -rf junk
	 mkdir junk && cd junk && ar x ../$(RAWDatabaseALIB) && ar x ../$(MDCALIB) && ar x ../$(ESDALIB) && ar x ../$(STEERBaseALIB) && ar r ../$(LIBPATH)/libAliMDC.a *.o && cd .. && rm -rf junk
	 $(LD) $(LDFLAGS) -o $(BINPATH)/alimdca $(alimdcCXXO) \
	 $(LIBPATH)/libAliMDC.a \
	 $(ROOTALIBDIR)/libRoot.a \
	 $(ROOTALIBDIR)/libfreetype.a $(ROOTALIBDIR)/libpcre.a \
	 -pthread -ldl -lcurses

alilibs-static: $(LIBPATH) modules $(ALLALIBS)

include  build/MakefileDA

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
	@echo "***** Cleaning up module.mk and temporary compilation files *****"
endif
	$(MUTE)rm -rf $(patsubst %,%/module.mk,$(MODULES))
	$(MUTE)rm -rf $(patsubst %,%/tgt_$(ALICE_TARGET),$(MODULES))

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
