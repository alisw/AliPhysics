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

FOBJS		= $(FSRCS:.F=.o)

# C Sources

CSRCS		= dummies.c

# C Objects

COBJS		= $(CSRCS:.c=.o) 

# C++ Sources

SRCS		= aliroot.cxx

# C++ Objects

OBJS		= $(SRCS:.cxx=.o) $(DICTO)

# Global Macros

ALLSRCS		= $(SRCS) $(FSRCS) $(CSRCS)
ALLOBJS		= $(OBJS) $(FOBJS) $(COBJS)

# C++ compilation flags

CXXFLAGS	= $(CXXOPTS) -I$(ROOTSYS)/include -I. -I$(ALICE_ROOT)/include
 
# FORTRAN compilation flags

FFLAGS		= $(FOPT)

##### TARGETS #####
 
# Target

PROGRAM		= $(ALICE_ROOT)/bin/aliroot

default:	$(PROGRAM)

depend:		$(ALLSRCS)

$(ALICE_ROOT)/bin/aliroot:	$(OBJS) $(COBJS) $(FOBJS)

LIBS = $(ALILIBS) $(G3LIBS) $(CLIBS) $(GLIBS) -lEG -lEGPythia -lpythia -lpdf

TOCLEAN		= $(PROGRAM) $(OBJS) $(COBJS) $(FOBJS)

############################### General Macros ################################

include $(ALICE_ROOT)/conf/GeneralMacros

############################ Dependencies #####################################

include Make-depend 
