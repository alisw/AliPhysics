############################### THijing Makefile ##############################

# Include machine specific definitions

include $(ALICE_ROOT)/conf/GeneralDef
include $(ALICE_ROOT)/conf/MachineDef.$(ALICE_TARGET)

PACKAGE = TIsajet

# FORTRAN sources

FSRCS	     = 

# C++ sources

SRCS          = TIsajet.cxx 

# C++ Headers

HDRS          = $(SRCS:.cxx=.h) TIsajetLinkDef.h

# Library dictionary

DICT          = TIsajetCint.cxx
DICTH         = $(DICT:.cxx=.h)
DICTO         = $(patsubst %.cxx,tgt_$(ALICE_TARGET)/%.o,$(DICT))

# Objects

OBJS          = $(patsubst %.cxx,tgt_$(ALICE_TARGET)/%.o,$(SRCS)) $(FOBJS) $(DICTO)

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

SLIBRARY	= $(LIBDIR)/libTIsajet.$(SL)
ALIBRARY	= $(LIBDIR)/libTIsajet.a

default:	$(SLIBRARY)

$(LIBDIR)/libTIsajet.$(SL):		$(OBJS)

$(DICT):				$(HDRS)

depend:			$(SRCS) $(FSRCS)

HEADERS		= TIsajet.h

TOCLEAN		= $(OBJS) *Cint.cxx *Cint.h

CHECKS        = $(patsubst %.cxx,check/%.viol,$(SRCS))

############################### General Macros ################################

include $(ALICE_ROOT)/conf/GeneralMacros

############################ Dependencies #####################################

-include tgt_$(ALICE_TARGET)/Make-depend 








