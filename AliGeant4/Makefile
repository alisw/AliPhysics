# $Id$

############################### AliGeant4 Makefile #############################

# Include Geant4 specific definitions

include $(G4INSTALL)/config/architecture.gmk

# Include machine specific definitions

include $(ALICE_ROOT)/conf/GeneralDef
include $(ALICE_ROOT)/conf/MachineDef.$(ALICE_TARGET)

# Include Alice Geant4 specific definitions

NEED_G4_LIBS := 1
include $(ALICE_ROOT)/conf/Geant4.gmk

PACKAGE = AliGeant4
TARGET  = aligeant4

# C++ sources

SRCS           = $(wildcard *.cxx)
SRCS          := $(filter-out $(TARGET).cxx,$(SRCS))

# C++ Headers

HDRS          =  $(wildcard *.icc) $(wildcard *.h)

# Library dictionary
# only for selected classes

DICT          = AliGeant4Cint.cxx
DICTH         = AliRunConfiguration.h AliGeant4LinkDef.h
DICTO         = $(patsubst %.cxx,tgt_$(ALICE_TARGET)/%.o,$(DICT))

# C++ Objects

OBJS          = $(patsubst %.cxx,tgt_$(ALICE_TARGET)/%.o,$(SRCS)) $(DICTO)

# Program C++ source, objects

PSRC          = $(TARGET).cxx
POBJ          = $(patsubst %.cxx,tgt_$(ALICE_TARGET)/%.o,$(PSRC))

# C++ compilation flags

CXXFLAGS      = $(CXXOPTS) -I$(ROOTSYS)/include -I. -I$(ALICE_ROOT)/include \
	        -I$(ALICE_ROOT)/TGeant4       

##### TARGETS #####
 
# Target

SLIBRARY	= $(LIBDIR)/libAliGeant4.$(SL)
ALIBRARY	= $(LIBDIR)/lib/libAliGeant4.a
PROGRAM		= $(BINDIR)/$(TARGET)

default:	$(SLIBRARY) $(PROGRAM)

$(LIBDIR)/libAliGeant4.$(SL):		$(OBJS)

$(BINDIR)/$(TARGET):			$(POBJ)

$(DICT):				$(DICTH)

depend:					$(SRCS)

#HEADERS =	$(filter-out %LinkDef.h,$(HDRS))

# filter-out Geant3 dependent libraries (ZDC requires symbols from G3)
# and add Geant4 dependent ones
ALILIBS := $(filter-out -lTGeant3,$(ALILIBS))
ALILIBS := $(filter-out -lZDC,$(ALILIBS))
ALILIBS += -lTGeant4 -lAliGeant4

LIBS = $(ALILIBS) $(G4LDLIBS) $(CLIBS) $(GLIBS) -lEG -lEGPythia -lpythia -lpdf

TOCLEAN	      = $(OBJS) *Cint.cxx *Cint.h

CHECKS        = $(patsubst %.cxx,check/%.viol,$(SRCS))

############################### General Macros ################################

include $(ALICE_ROOT)/conf/GeneralMacros


# Geant4 specific compiler flags

CXXFLAGS += $(CPPFLAGS)


############################ Dependencies #####################################

-include tgt_$(ALICE_TARGET)/Make-depend 
