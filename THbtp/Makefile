############################### THBTprocessor Makefile ##############################

# Include machine specific definitions

include $(ALICE_ROOT)/conf/GeneralDef
include $(ALICE_ROOT)/conf/MachineDef.$(ALICE_TARGET)




PACKAGE = THbtp

LIBS = $(ALILIBS) -lHBTP

# FORTRAN sources

FSRCS	     = 

# C++ sources

SRCS          = THBTprocessor.cxx 

# C++ Headers

HDRS          = $(SRCS:.cxx=.h) THBTprocessorLinkDef.h 

# Library dictionary

DICT          = THBTprocessorCint.cxx
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

SLIBRARY	= $(LIBDIR)/libTHbtp.$(SL)
ALIBRARY	= $(LIBDIR)/libTHbtp.a

default:	$(SLIBRARY)

$(LIBDIR)/libTHbtp.$(SL):		$(OBJS)

$(DICT):				$(HDRS)

depend:			$(SRCS) $(FSRCS)

HEADERS		= HBTprocCOMMON.h THBTprocessor.h 

TOCLEAN		= $(OBJS) *Cint.cxx *Cint.h

CHECKS        = $(patsubst %.cxx,check/%.viol,$(SRCS))

############################### General Macros ################################

include $(ALICE_ROOT)/conf/GeneralMacros

############################ Dependencies #####################################

-include tgt_$(ALICE_TARGET)/Make-depend 







