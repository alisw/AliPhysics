############################### RICH Makefile #################################

# Include machine specific definitions

include $(ALICE_ROOT)/conf/GeneralDef
include $(ALICE_ROOT)/conf/MachineDef.$(ALICE_TARGET)

PACKAGE = RICH

# C++ sources

SRCS          = AliRICH.cxx

# C++ Headers

HDRS          = $(SRCS:.cxx=.h) RICHLinkDef.h

# Library dictionary

DICT          = RICHCint.cxx
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

CXXFLAGS      += -I$(ALICE_ROOT)/TGeant3
DEPINC        += -I$(ALICE_ROOT)/TGeant3
 
# FORTRAN compilation flags

FFLAGS      = $(FOPT)

##### TARGETS #####
 
# Target

SLIBRARY	= $(LIBDIR)/libRICH.$(SL)
ALIBRARY	= $(LIBDIR)/libRICH.a

default:	$(SLIBRARY)

$(LIBDIR)/libRICH.$(SL):		$(OBJS)

$(DICT):				$(HDRS)

depend:					$(SRCS)

TOCLEAN		= $(OBJS) *Cint.cxx *Cint.h

############################### General Macros ################################

include $(ALICE_ROOT)/conf/GeneralMacros

############################ Dependencies #####################################

-include Make-depend 
