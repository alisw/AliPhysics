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
DICTO         = $(patsubst %.cxx,tgt_$(ALICE_TARGET)/%.o,$(DICT))

# Fortran sources

FSRCS	      = algo.F

# FORTRAN Objectrs

FOBJS         = $(patsubst %.F,tgt_$(ALICE_TARGET)/%.o,$(FSRCS))

# C Objects

COBJS         = $(patsubst %.c,tgt_$(ALICE_TARGET)/%.o,$(CSRCS))

# C++ Objects

OBJS          = $(patsubst %.cxx,tgt_$(ALICE_TARGET)/%.o,$(SRCS)) $(DICTO)

# C++ compilation flags

CXXFLAGS      = $(CXXOPTS) -I$(ROOTSYS)/include -I. -I$(ALICE_ROOT)/include/
 
# FORTRAN compilation flags

FFLAGS      = $(FOPT)

##### TARGETS #####
 
# Target

SLIBRARY	= $(LIBDIR)/libMUON.$(SL)
ALIBRARY	= $(LIBDIR)/libMUON.a

default:	$(SLIBRARY)

$(LIBDIR)/libMUON.$(SL):	$(OBJS) $(FOBJS)

$(DICT):			$(HDRS)

depend:				$(SRCS)

TOCLEAN		= $(OBJS) $(FOBJS) *Cint.cxx *Cint.h

############################### General Macros ################################

include $(ALICE_ROOT)/conf/GeneralMacros

############################ Dependencies #####################################

-include tgt_$(ALICE_TARGET)/Make-depend 
 
