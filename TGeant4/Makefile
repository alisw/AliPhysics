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
DICTH         = TGeant4.h TG4VRunConfiguration.h TGeant4LinkDef.h
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


############################ Dependencies #####################################

-include tgt_$(ALICE_TARGET)/Make-depend 
