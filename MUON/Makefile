############################### MUON Makefile #################################

# Include machine specific definitions

include $(ALICE_ROOT)/conf/GeneralDef
include $(ALICE_ROOT)/conf/MachineDef.$(ALICE_TARGET)

PACKAGE = MUON

# C++ sources

SRCS          = AliMUON.cxx AliMUONdisplay.cxx AliMUONpoints.cxx \
                AliMUONsegmentv1.cxx AliMUONv0.cxx AliMUONv01.cxx

# C++ Headers

HDRS          = $(SRCS:.cxx=.h) AliMUONConst.h MUONLinkDef.h

# Library dictionary

DICT          = MUONCint.cxx
DICTH         = $(DICT:.cxx=.h)
DICTO         = $(DICT:.cxx=.o)

# Fortran sources

FSRCS	      = algo.F

# FORTRAN Objectrs

FOBJS         = $(FSRCS:.F=.o)

# C Objects

COBJS         = $(CSRCS:.c=.o)

# C++ Objects

OBJS          = $(SRCS:.cxx=.o) $(DICTO)

# C++ compilation flags

CXXFLAGS      = $(CXXOPTS) -I$(ROOTSYS)/include -I. -I$(ALICE_ROOT)/include/
 
# FORTRAN compilation flags

FFLAGS      = $(FOPT)

##### TARGETS #####
 
# Target

SLIBRARY	= $(ALICE_ROOT)/lib/libMUON.$(SL)
ALIBRARY	= $(ALICE_ROOT)/lib/libMUON.a

default:	$(SLIBRARY)

$(ALICE_ROOT)/lib/libMUON.$(SL):	$(OBJS) $(FOBJS)

$(DICT):			$(HDRS)

depend:				$(SRCS)

TOCLEAN		= $(OBJS) $(FOBJS) *Cint.cxx *Cint.h

############################### General Macros ################################

include $(ALICE_ROOT)/conf/GeneralMacros

############################ Dependencies #####################################

include Make-depend 
 
