############################### TRD Makefile ##################################

# Include machine specific definitions

include $(ALICE_ROOT)/conf/GeneralDef
include $(ALICE_ROOT)/conf/MachineDef.$(ALICE_TARGET)

# Package name

PACKAGE	      = TRD

# C++ sources

SRCS          = AliTRD.cxx	AliTRDv0.cxx      AliTRDv1.cxx    AliTRDv2.cxx 

# C++ Headers

HDRS          = $(SRCS:.cxx=.h) TRDLinkDef.h

# Library dictionary

DICT          = TRDCint.cxx
DICTH         = $(DICT:.cxx=.h)
DICTO         = $(DICT:.cxx=.o)

# FORTRAN Objectrs

FOBJS         = $(FSRCS:.f=.o)

# C Objects

COBJS         = $(CSRCS:.c=.o)

# C++ Objects

OBJS          = $(SRCS:.cxx=.o) $(DICTO)

# C++ compilation flags

CXXFLAGS      = $(CXXOPTS) -I$(ROOTSYS)/include -I. -I$(ALICE_ROOT)/include
 
# FORTRAN compilation flags

FFLAGS      = $(FOPT)

##### TARGETS #####
 
# Target

SLIBRARY	= $(ALICE_ROOT)/lib/libTRD.$(SL)
ALIBRARY	= $(ALICE_ROOT)/lib/libTRD.a

default:	$(SLIBRARY)

$(ALICE_ROOT)/lib/libTRD.$(SL):	$(OBJS)

TRDCint.cxx:			$(HDRS)

depend:				$(SRCS)

TOCLEAN			= $(OBJS) *Cint.h *Cint.cxx

############################### General Macros ################################

include $(ALICE_ROOT)/conf/GeneralMacros

############################ Dependencies #####################################

-include Make-depend 
