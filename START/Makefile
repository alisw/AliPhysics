
###############################START Makefile##################################

# Include machine specific definitions

include $(ALICE_ROOT)/conf/GeneralDef
include $(ALICE_ROOT)/conf/MachineDef.$(ALICE_TARGET)

PACKAGE = START

# C++ sources

SRCS          = AliSTART.cxx AliSTARTv0.cxx AliSTARTv1.cxx AliSTARThit.cxx AliSTARTdigit.cxx AliSTARTvertex.cxx

# C++ Headers

HDRS          = $(SRCS:.cxx=.h) STARTLinkDef.h

# Library dictionary

DICT          = STARTCint.cxx
DICTH         = $(DICT:.cxx=.h)
DICTO         = $(patsubst %.cxx,tgt_$(ALICE_TARGET)/%.o,$(DICT))

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

SLIBRARY	= $(LIBDIR)/libSTART.$(SL)
ALIBRARY	= $(LIBDIR)/libSTART.a

default:	$(SLIBRARY)

$(LIBDIR)/libSTART.$(SL):	$(OBJS)

$(DICT):			$(HDRS)

depend:				$(SRCS)

TOCLEAN		 = $(OBJS) *Cint.cxx *Cint.h

CHECKS        = $(patsubst %.cxx,check/%.viol,$(SRCS))

############################### General Macros ################################

include $(ALICE_ROOT)/conf/GeneralMacros

############################ Dependencies #####################################

-include tgt_$(ALICE_TARGET)/Make-depend 

