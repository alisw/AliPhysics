# $Id$

############################### TGeant4 Makefile ###############################

# Include Geant4 specific definitions

include $(G4INSTALL)/config/architecture.gmk

# Include machine specific definitions

include $(ALICE_ROOT)/conf/GeneralDef
include $(ALICE_ROOT)/conf/MachineDef.$(ALICE_TARGET)

# Include Alice Geant4 specific definitions

include $(ALICE_ROOT)/conf/Geant4.gmk

PACKAGE = TGeant4

# C++ sources

SRCS           = $(wildcard *.cxx)

# C++ Headers

HDRS          =  $(wildcard *.icc) $(wildcard *.h)

# Library Root dictionary
# only for selected classes
# (temporarily excluded for OSF1 due to compilation errors)

ifneq ($(ALICE_TARGET),OSF1)
DICT          = TGeant4Cint.cxx
DICTH         = TGeant4.h TG4VRunConfiguration.h \
	        TG4Editor.h TG4GUI.h TG4GeometryGUI.h TG4GuiVolume.h \
	        TGeant4LinkDef.h
DICTO         = $(patsubst %.cxx,tgt_$(ALICE_TARGET)/%.o,$(DICT))
endif

# C++ Objects

OBJS          = $(patsubst %.cxx,tgt_$(ALICE_TARGET)/%.o,$(SRCS)) $(DICTO)

# C++ compilation flags

CXXFLAGS      = $(CXXOPTS) -I$(ROOTSYS)/include -I. -I$(ALICE_ROOT)/include/
 
 
##### TARGETS #####
 
# Target

SLIBRARY	= $(LIBDIR)/libTGeant4.$(SL)
ALIBRARY	= $(LIBDIR)/libTGeant4.a

default:	$(SLIBRARY)

$(LIBDIR)/libTGeant4.$(SL):		$(OBJS)
$(LIBDIR)/libTGeant4.a:			$(OBJS)

$(DICT):	                        $(DICTH)

depend:					$(SRCS)

clean: macroclean

#HEADERS =	$(filter-out %LinkDef.h,$(HDRS))

TOCLEAN	=       $(OBJS) *Cint.cxx *Cint.h

CHECKS  =       $(patsubst %.cxx,check/%.viol,$(SRCS))

############################### General Macros ################################

include $(ALICE_ROOT)/conf/GeneralMacros

# Geant4 specific compiler flags

CXXFLAGS += $(CPPFLAGS)

# Remove optimisation for TG4PhysicsConstructorHadron

ifeq ($(ALICE_TARGET),Linux)
  COMPILER := $(shell g++ -v 2>&1 | grep version | grep egcs)
  ifeq ($(findstring 1.1,$(COMPILER)),1.1)
    CXXFLAGS_WITHOUT_O := $(filter-out -O% -Wuninitialized, $(CXXFLAGS))
    tgt_$(ALICE_TARGET)/TG4PhysicsConstructorHadron.o : TG4PhysicsConstructorHadron.cxx
      ifeq ($(CXXVERBOSE),)
	@$(CXX) $(CXXFLAGS_WITHOUT_O) -c -o $@ TG4PhysicsConstructorHadron.cxx
      else
	$(CXX) $(CXXFLAGS_WITHOUT_O) -c -o $@ TG4PhysicsConstructorHadron.cxx
      endif
  endif
endif      

ifeq ($(ALICE_TARGET),HP-UX.aCC)
  # remove optimisation for TG4PhysicsConstructorHadron
  # to be tested
  COMPILER := $(shell aCC -V 2>&1)
  ifeq ($(COMPILER), aCC: HP ANSI C++ B3910B A.01.15)
    CXXFLAGS_WITHOUT_O := $(filter-out -O% , $(CXXFLAGS))
    CXXFLAGS_WITHOUT_O := $(filter-out +O% , $(CXXFLAGS_WITHOUT_O))  
    tgt_$(ALICE_TARGET)/TG4PhysicsConstructorHadron.o : TG4PhysicsConstructorHadron.cxx
      ifeq ($(CXXVERBOSE),)
	@$(CXX) $(CXXFLAGS_WITHOUT_O) -c $*.cxx -o tgt_$(ALICE_TARGET)/$*.o
      else
	$(CXX) $(CXXFLAGS_WITHOUT_O) -c $*.cxx -o tgt_$(ALICE_TARGET)/$*.o
      endif
  endif    
endif      


############################ Dependencies #####################################

-include tgt_$(ALICE_TARGET)/Make-depend 
