############################### ZDC Makefile ##################################

# Include machine specific definitions

include $(ALICE_ROOT)/conf/GeneralDef
include $(ALICE_ROOT)/conf/MachineDef.$(ALICE_TARGET)

PACKAGE = ZDC

# FORTRAN sources

FSRCS	      = AliZDCf.F

# C++ sources

SRCS          = AliZDC.cxx

# Fortran sources

FSRCS	       = AliZDCf.F 

# C++ Headers

HDRS          = $(SRCS:.cxx=.h) ZDCLinkDef.h

# Library dictionary

DICT          = ZDCCint.cxx
DICTH         = $(DICT:.cxx=.h)
DICTO         = $(patsubst %.cxx,$(ALICE_TARGET)/%.o,$(DICT))

# FORTRAN Objectrs

FOBJS         = $(patsubst %.F,$(ALICE_TARGET)/%.o,$(FSRCS))

# C Objects

COBJS         = $(patsubst %.c,$(ALICE_TARGET)/%.o,$(CSRCS)

# C++ Objects

OBJS          = $(patsubst %.cxx,$(ALICE_TARGET)/%.o,$(SRCS)) $(DICTO)

# Global Macros

ALLOBJS       = $(OBJS) $(FOBJS)

# C++ compilation flags

CXXFLAGS      = $(CXXOPTS) -I$(ROOTSYS)/include -I. -I$(ALICE_ROOT)/include
 
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

############################### General Macros ################################

include $(ALICE_ROOT)/conf/GeneralMacros

############################ Dependencies #####################################

-include Make-depend 

