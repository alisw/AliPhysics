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
DICTO         = $(DICT:.cxx=.o)

# FORTRAN Objectrs

FOBJS         = $(FSRCS:.F=.o)

# C Objects

COBJS         = $(CSRCS:.c=.o)

# C++ Objects

OBJS          = $(SRCS:.cxx=.o) $(DICTO)

# Global Macros

ALLOBJS       = $(OBJS) $(FOBJS)

# C++ compilation flags

CXXFLAGS      = $(CXXOPTS) -I$(ROOTSYS)/include -I. -I$(ALICE_ROOT)/include
 
# FORTRAN compilation flags

FFLAGS      = $(FOPT) -I$(ALICE_ROOT)/GEANT321

##### TARGETS #####
 
# Target

SLIBRARY	= $(ALICE_ROOT)/lib/libZDC.$(SL)
ALIBRARY	= $(ALICE_ROOT)/lib/libZDC.a

default:	$(SLIBRARY)

$(ALICE_ROOT)/lib/libZDC.$(SL):	$(OBJS) $(FOBJS)

$(DICT):			$(HDRS)

DEPINC	+=	-I$(ALICE_ROOT)/GEANT321

depend:				$(SRCS) $(FSRCS)

TOCLEAN		= $(OBJS) $(FOBJS) *Cint.cxx *Cint.h

############################### General Macros ################################

include $(ALICE_ROOT)/conf/GeneralMacros

############################ Dependencies #####################################

include Make-depend 

