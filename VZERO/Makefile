###############################VZERO Makefile##################################

# Include machine specific definitions

include $(ALICE_ROOT)/conf/GeneralDef
include $(ALICE_ROOT)/conf/MachineDef.$(ALICE_TARGET)

PACKAGE = VZERO

# C++ sources

SRCS          = AliVZERO.cxx AliVZEROv0.cxx AliVZEROhit.cxx AliVZEROdigit.cxx

# C++ Headers

HDRS          = $(SRCS:.cxx=.h) VZEROLinkDef.h

# Library dictionary

DICT          = VZEROCint.cxx
DICTH         = $(DICT:.cxx=.h)
DICTO         = $(patsubst %.cxx,tgt_$(ALICE_TARGET)/%.o,$(DICT))

# FORTRAN Objectrs

FOBJS         = $(patsubst %.F,tgt_$(ALICE_TARGET)/%.o,$(FSRCS))

# C Objects

COBJS         = $(patsubst %.c,tgt_$(ALICE_TARGET)/%.o,$(CSRCS))

# C++ Objects

OBJS          = $(patsubst %.cxx,tgt_$(ALICE_TARGET)/%.o,$(SRCS)) $(DICTO)

# C++ compilation flags

CXXFLAGS      = $(CXXOPTS) -I$(ROOTSYS)/include -I. -I$(ALICE_ROOT)/include

CXXFLAGS      += -I$(ALICE_ROOT)/TGeant3 -I$(ALICE_ROOT)/STRUCT
DEPINC        += -I$(ALICE_ROOT)/TGeant3
 
# FORTRAN compilation flags

FFLAGS      = $(FOPT)

##### TARGETS #####
 
# Target

SLIBRARY	= $(LIBDIR)/libVZERO.$(SL)
ALIBRARY	= $(LIBDIR)/libVZERO.a

default:	$(SLIBRARY)

$(LIBDIR)/libVZERO.$(SL):	$(OBJS)

$(DICT):			$(HDRS)

depend:				$(SRCS)

TOCLEAN		 = $(OBJS) *Cint.cxx *Cint.h

############################### General Macros ################################

include $(ALICE_ROOT)/conf/GeneralMacros

############################ Dependencies #####################################

-include tgt_$(ALICE_TARGET)/Make-depend 
 

### Target check creates violation reports (.viol), which depend on
### stripped files (.ii), which in turn depend on preprocessed
### files (.i). Dependences are in conf/GeneralDef.

CHECKS        = $(patsubst %.cxx,check/%.viol,$(SRCS))

check:          $(CHECKS)

# DO NOT DELETE
