############################### CASTOR Makefile ###############################

# Include machine specific definitions

include $(ALICE_ROOT)/conf/GeneralDef
include $(ALICE_ROOT)/conf/MachineDef.$(ALICE_TARGET)

PACKAGE = CASTOR

# C++ sources

SRCS          = AliCASTOR.cxx

# C++ Headers

HDRS          = $(SRCS:.cxx=.h) CASTORLinkDef.h

# Library dictionary

DICT          = CASTORCint.cxx
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

SLIBRARY	= $(ALICE_ROOT)/lib/libCASTOR.$(SL)
ALIBRARY	= $(ALICE_ROOT)/lib/libCASTOR.a

default:	$(SLIBRARY)

$(ALICE_ROOT)/lib/libCASTOR.$(SL):		$(OBJS)

$(DICT):				$(HDRS)

depend:					$(SRCS)

TOCLEAN		= $(OBJS) *Cint.cxx *Cint.h

############################### General Macros ################################

include $(ALICE_ROOT)/conf/GeneralMacros

############################ Dependencies #####################################

-include Make-depend 
