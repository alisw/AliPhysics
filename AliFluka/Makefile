# $Id$

############################### AliFluka Makefile #############################

# Include Geant4 specific definitions

include $(FLUGGINSTALL)/config/architecture.gmk

# Include machine specific definitions

include $(ALICE_ROOT)/conf/GeneralDef
include $(ALICE_ROOT)/conf/MachineDef.$(ALICE_TARGET)

# Include Alice Geant4 specific definitions

NEED_G4_LIBS := 1
include $(ALICE_ROOT)/conf/Fluka.gmk

PACKAGE = AliFluka
TARGET  = alifluka

# C++ sources

SRCS           = $(wildcard *.cxx)
SRCS          := $(filter-out $(TARGET).cxx,$(SRCS))

# C++ Headers

HDRS          =  $(wildcard *.icc) $(wildcard *.h) $(wildcard *.hh)

# Library Root dictionary
# only for selected classes
# (temporarily excluded for OSF1 due to compilation errors)

ifneq ($(ALICE_TARGET),OSF1)
#DICT          = AliFlukaCint.cxx
#DICTH         = AliFlukaLinkDef.h
#DICTO         = $(patsubst %.cxx,tgt_$(ALICE_TARGET)/%.o,$(DICT))
endif

# C++ Objects

OBJS          = $(patsubst %.cxx,tgt_$(ALICE_TARGET)/%.o,$(SRCS)) $(DICTO)

# Program C++ source, objects

PSRC          = $(TARGET).cxx
POBJ          = $(patsubst %.cxx,tgt_$(ALICE_TARGET)/%.o,$(PSRC))

# C++ compilation flags

CXXFLAGS      = $(CXXOPTS) -I$(ROOTSYS)/include -I. -I$(ALICE_ROOT)/include \
	        -I$(ALICE_ROOT)/TGeant4 -I$(ALICE_ROOT)/AliGeant4 \
                -I$(ALICE_ROOT)/STRUCT   

##### TARGETS #####
 
# Target

SLIBRARY	= $(LIBDIR)/libAliFluka.$(SL)
ALIBRARY	= $(LIBDIR)/libAliFluka.a
PROGRAM		= $(BINDIR)/$(TARGET)

default:	$(SLIBRARY) $(PROGRAM)

$(LIBDIR)/libAliFluka.$(SL):		$(OBJS)
$(LIBDIR)/libAliFluka.a:		$(OBJS)

$(BINDIR)/$(TARGET):			$(POBJ)

$(DICT):				$(DICTH)

depend:					$(SRCS)

clean: macroclean

#HEADERS =	$(filter-out %LinkDef.h,$(HDRS))

# filter-out Geant3 dependent libraries (ZDC requires symbols from G3)
# and add Geant4 dependent ones
ALILIBS := $(filter-out -lTGeant3,$(ALILIBS))
ALILIBS := $(filter-out -lZDC,$(ALILIBS))
ALILIBS := $(filter-out -lITS,$(ALILIBS))
ALILIBS += -lAliFluka -lAliGeant4_1 -lTGeant4_1

LIBS = $(ALILIBS) $(CLIBS) $(GLIBS) -lEG -lEGPythia6 -lpythia6 -lpdf \
	-lTHijing -lhijing  -lTMevSim -lmevsim   -lTHbtp -lHBTP \
	$(G4LDLIBS)

TOCLEAN	      = $(OBJS) $(POBJ) *Cint.cxx *Cint.h

CHECKS        = $(patsubst %.cxx,check/%.viol,$(SRCS))

############################### General Macros ################################

include $(ALICE_ROOT)/conf/GeneralMacros


# Geant4 specific compiler flags

CXXFLAGS += $(CPPFLAGS)


############################ Dependencies #####################################

-include tgt_$(ALICE_TARGET)/Make-depend 
