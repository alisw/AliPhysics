############################### ALIROOT Makefile ##############################

# Include machine specific definitions

include $(ALICE_ROOT)/conf/GeneralDef
include $(ALICE_ROOT)/conf/MachineDef.$(ALICE_TARGET)

PACKAGE		= ALIROOT

# C++ sources

SRCS		= 

# C++ Headers

HDRS		=

# Library dictionary

DICT		= 
DICTH		= $(DICT:.cxx=.h)
DICTO		= $(DICT:.cxx=.o)

# FORTRAN Sources

FSRCS		= loadmod.F pydata.F ludata.F

# FORTRAN Objectrs

FOBJS		= $(patsubst %.F,tgt_$(ALICE_TARGET)/%.o,$(FSRCS))

# C Sources

CSRCS		= dummies.c

# C Objects

COBJS		= $(patsubst %.c,tgt_$(ALICE_TARGET)/%.o,$(CSRCS)) 

# C++ Sources

SRCS		= aliroot.cxx

# C++ Objects

OBJS		= $(patsubst %.cxx,tgt_$(ALICE_TARGET)/%.o,$(SRCS)) $(DICTO)

# Global Macros

ALLSRCS		= $(SRCS) $(FSRCS) $(CSRCS)
ALLOBJS		= $(OBJS) $(FOBJS) $(COBJS)

# C++ compilation flags

CXXFLAGS	= $(CXXOPTS) -I$(ROOTSYS)/include -I. -I$(ALICE_ROOT)/include
 
# FORTRAN compilation flags

FFLAGS		= $(FOPT)

##### TARGETS #####
 
# Target

PROGRAM		= $(BINDIR)/aliroot

default:	$(PROGRAM)

depend:		$(ALLSRCS)

$(BINDIR)/aliroot:	$(OBJS) $(COBJS) $(FOBJS)

LIBS = $(ALILIBS) $(G3LIBS) $(CLIBS) $(GLIBS) -lEG -lEGPythia -lpythia -lpdf

TOCLEAN		= $(PROGRAM) $(OBJS) $(COBJS) $(FOBJS)

############################### General Macros ################################

include $(ALICE_ROOT)/conf/GeneralMacros

############################ Dependencies #####################################

-include tgt_$(ALICE_TARGET)/Make-depend 















### Target check creates violation reports (.viol), which depend on
### stripped files (.ii), which in turn depend on preprocessed
### files (.i). Dependences are in conf/GeneralDef.

CHECKS        = $(patsubst %.cxx,check/%.viol,$(SRCS))

check:          $(CHECKS)







