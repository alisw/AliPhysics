############################### TPC Makefile ##################################

# Include machine specific definitions

include $(ALICE_ROOT)/conf/GeneralDef
include $(ALICE_ROOT)/conf/MachineDef.$(ALICE_TARGET)

PACKAGE = TPC

# C++ sources

SRCS          = AliTPC.cxx  AliTPCv0.cxx    AliTPCv1.cxx      AliTPCv2.cxx  \
		AliTPCD.cxx AliTPCPRF2D.cxx AliTPCRF1D.cxx AliTPCParam.cxx

# C++ Headers

HDRS          = $(SRCS:.cxx=.h) AliTPCSecGeo.h TPCLinkDef.h

# Library dictionary

DICT          = TPCCint.cxx
DICTH         = $(DICT:.cxx=.h)
DICTO         = $(DICT:.cxx=.o)

# FORTRAN Objectrs

FOBJS         = $(FSRCS:.f=.o)

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

SLIBRARY	= $(ALICE_ROOT)/lib/libTPC.$(SL)

default:	$(SLIBRARY)

$(ALICE_ROOT)/lib/libTPC.$(SL):  $(OBJS)

$(DICT):		$(HDRS)

depend:			$(SRCS)

TOCLEAN			= $(OBJS) *Cint.h *Cint.cxx

############################### General Macros ################################

include $(ALICE_ROOT)/conf/GeneralMacros

############################ Dependencies #####################################

-include Make-depend 
