############################### CRT Makefile ##################################

# Include machine specific definitions

include $(ALICE_ROOT)/conf/GeneralDef
include $(ALICE_ROOT)/conf/MachineDef.$(ALICE_TARGET)

PACKAGE = CRT

# C++ sources

SRCS  = AliCRT.cxx AliCRTv0.cxx AliCRTv1.cxx AliGenCRT.cxx AliCRThit.cxx AliCRTdigit.cxx AliCRTConstants.cxx 

# C++ Headers

HDRS          = $(SRCS:.cxx=.h) CRTLinkDef.h

# Library dictionary

DICT          = CRTCint.cxx
DICTH         = $(DICT:.cxx=.h)
DICTO         = $(patsubst %.cxx,tgt_$(ALICE_TARGET)/%.o,$(DICT))

# FORTRAN Objectrs

FOBJS         = $(patsubst %.f,tgt_$(ALICE_TARGET)/%.o,$(FSRCS))

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

SLIBRARY	= $(LIBDIR)/libCRT.$(SL)
ALIBRARY	= $(LIBDIR)/libCRT.a

default:	$(SLIBRARY)

$(LIBDIR)/libCRT.$(SL):	$(OBJS)

$(DICT):			$(HDRS)

depend:				$(SRCS)

TOCLEAN			= $(OBJS) *Cint.cxx *Cint.h

CHECKS		        = $(patsubst %.cxx,check/%.viol,$(SRCS))

############################### General Macros ################################

include $(ALICE_ROOT)/conf/GeneralMacros

############################ Dependencies #####################################

-include tgt_$(ALICE_TARGET)/Make-depend 
# DO NOT DELETE
