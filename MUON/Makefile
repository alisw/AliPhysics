############################### MUON Makefile #################################

# Include machine specific definitions

include $(ALICE_ROOT)/conf/GeneralDef
include $(ALICE_ROOT)/conf/MachineDef.$(ALICE_TARGET)

PACKAGE = MUON

# C++ sources

SRCS          = AliMUONchamber.cxx AliMUONSegResV0.cxx AliMUONSegResV01.cxx \
               AliMUONSegResV02.cxx AliMUONSegResV04.cxx AliMUONSegResV05.cxx\
               AliMUONSegResV1.cxx AliMUON.cxx AliMUONv0.cxx \
               AliMUONdisplay.cxx AliMUONpoints.cxx \
               AliMUONClusterFinder.cxx  AliMUONClusterFinderv0.cxx \
               AliMUONHitMap.cxx AliMUONTUBE.cxx

# C++ Headers

HDRS          = $(SRCS:.cxx=.h) AliMUONConst.h MUONLinkDef.h

# Library dictionary

DICT          = MUONCint.cxx
DICTH         = $(DICT:.cxx=.h)
DICTO         = $(patsubst %.cxx,tgt_$(ALICE_TARGET)/%.o,$(DICT))

# Fortran sources

FSRCS	      = algo.F reco_muon.F

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
 

### Target check creates violation reports (.viol), which depend on
### stripped files (.ii), which in turn depend on preprocessed
### files (.i). Dependences are in conf/GeneralDef.

CHECKS        = $(patsubst %.cxx,check/%.viol,$(SRCS))

check:          $(CHECKS)

