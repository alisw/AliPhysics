############################### ZDC Makefile ##################################

# Include machine specific definitions

include $(ALICE_ROOT)/conf/GeneralDef
include $(ALICE_ROOT)/conf/MachineDef.$(ALICE_TARGET)

PACKAGE = ZDC

# FORTRAN sources

FSRCS	      =

# C++ sources

SRCS          = AliZDC.cxx AliZDCv1.cxx AliZDCHit.cxx AliZDCDigit.cxx AliGenZDC.cxx

# C++ Headers

HDRS          = $(SRCS:.cxx=.h) $(ROOTSYS)/include/TClonesArray.h ZDCLinkDef.h

# Library dictionary

DICT          = ZDCCint.cxx
DICTH         = $(DICT:.cxx=.h)
DICTO         = $(patsubst %.cxx,tgt_$(ALICE_TARGET)/%.o,$(DICT))

# FORTRAN Objectrs

FOBJS         = $(patsubst %.F,tgt_$(ALICE_TARGET)/%.o,$(FSRCS))

# C Objects

COBJS         = $(patsubst %.c,tgt_$(ALICE_TARGET)/%.o,$(CSRCS)

# C++ Objects

OBJS          = $(patsubst %.cxx,tgt_$(ALICE_TARGET)/%.o,$(SRCS)) $(DICTO)

# Global Macros

ALLOBJS       = $(OBJS) $(FOBJS)

# C++ compilation flags

CXXFLAGS      = $(CXXOPTS) -I$(ROOTSYS)/include -I. -I$(ALICE_ROOT)/include \
		-I$(ALICE_ROOT)/RALICE
 
# FORTRAN compilation flags

FFLAGS      = $(FOPT) -I$(ALICE_ROOT)/GEANT321

##### TARGETS #####
 
# Target

SLIBRARY	= $(LIBDIR)/libZDC.$(SL)
ALIBRARY	= $(LIBDIR)/libZDC.a

default:	$(SLIBRARY)

$(LIBDIR)/libZDC.$(SL):	$(OBJS) $(FOBJS)

$(DICT):			$(HDRS)

DEPINC	+=	-I$(ALICE_ROOT)/GEANT321

depend:				$(SRCS) $(FSRCS)

TOCLEAN		= $(OBJS) $(FOBJS) *Cint.cxx *Cint.h

CHECKS        = $(patsubst %.cxx,check/%.viol,$(SRCS))

############################### General Macros ################################

include $(ALICE_ROOT)/conf/GeneralMacros

############################ Dependencies #####################################

-include tgt_$(ALICE_TARGET)/Make-depend 







