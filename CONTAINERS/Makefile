# Include machine specific definitions

include $(ALICE_ROOT)/conf/GeneralDef
include $(ALICE_ROOT)/conf/MachineDef.$(ALICE_TARGET)

PACKAGE = CONTAINERS

# C++ sources 

SRCS          =  AliMemArray.cxx AliClassInfo.cxx  AliDataType.cxx\
		 AliObjectArray.cxx    TArrayOfArray.cxx  AliArrayBranch.cxx\
		AliSegmentID.cxx  AliSegmentArray.cxx \
		AliDigits.cxx AliArrayI.cxx AliArrayS.cxx
# C++ Headers

HDRS          = $(SRCS:.cxx=.h) CONTAINERSLinkDef.h

# Library dictionary

DICT          = CONTAINERSCint.cxx
DICTH         = $(DICT:.cxx=.h)
DICTO         = $(patsubst %.cxx,tgt_$(ALICE_TARGET)/%.o,$(DICT))

# FORTRAN Objectrs

FOBJS         = $(FSRCS:.f=.o)
# C Objects

COBJS         = $(CSRCS:.c=.o)

# C++ Objects

OBJS          = $(patsubst %.cxx,tgt_$(ALICE_TARGET)/%.o,$(SRCS)) $(DICTO)

# C++ compilation flags

CXXFLAGS      = $(CXXOPTS) -g -I$(ALICE_ROOT)/CONTAINERS  -I$(ROOTSYS)/include -I. -I$(ALICE_ROOT)/include -I$(ALICE_ROOT)/TPC 
#CXXFLAGS      = $(CXXOPTS) -I$(ROOTSYS)/include -I. -I$(ALICE_ROOT)/include/
 
# FORTRAN compilation flags

FFLAGS      = $(FOPT)
 
##### TARGETS #######
 
# Target

SLIBRARY	= $(LIBDIR)/libCONTAINERS.$(SL)

default:	$(SLIBRARY)

$(LIBDIR)/libCONTAINERS.$(SL):  $(OBJS)

$(DICT):		$(HDRS)

depend:			$(SRCS)

TOCLEAN			= $(OBJS) *Cint.h *Cint.cxx

CHECKS        = $(patsubst %.cxx,check/%.viol,$(SRCS))

############################### General Macros ################################

include $(ALICE_ROOT)/conf/GeneralMacros

############################ Dependencies #####################################

-include tgt_$(ALICE_TARGET)/Make-depend 


