############################### TGeant3 Makefile ##############################

# Include machine specific definitions

include $(ALICE_ROOT)/conf/GeneralDef
include $(ALICE_ROOT)/conf/MachineDef.$(ALICE_TARGET)

PACKAGE = TGeant3

# FORTRAN sources

FSRCS	     = galicef.F

# C++ sources

GSRCS          = TGeant3.cxx TPaveTree.cxx THIGZ.cxx gucode.cxx

# C++ Headers

GHDRS          = $(filter-out gucode.h,$(GSRCS:.cxx=.h)) TGeant3LinkDef.h

# Library dictionary

GDICT          = TGeant3Cint.cxx
GDICTH         = $(GDICT:.cxx=.h)
GDICTO         = $(GDICT:.cxx=.o)

# Geant Objects

GOBJS          = $(GSRCS:.cxx=.o) $(FOBJS) $(GDICTO)

# FORTRAN Objectrs

FOBJS         = $(FSRCS:.F=.o)

# C Objects

COBJS         = $(CSRCS:.c=.o)

# C++ Objects

OBJS          = $(SRCS:.cxx=.o) $(DICTO)

# Dummy Geant321 class

# C++ sources

DSRCS          = TGeant3Dummy.cxx

# C++ Headers

DHDRS          = TGeant3.h TGeant3LinkDef.h

# Library dictionary

DDICT          = TGeant3DummyCint.cxx
DDICTH         = $(DDICT:.cxx=.h)
DDICTO         = $(DDICT:.cxx=.o)

# Geant Objects

DOBJS          = $(DSRCS:.cxx=.o) $(DDICTO)

# Dummy Routines in Geant321

JSRCS		= G3Dummy.cxx
JOBJS		= G3Dummy.o

# C++ compilation flags

CXXFLAGS      = $(CXXOPTS) -I$(ROOTSYS)/include -I. -I$(ALICE_ROOT)/include/
 
# FORTRAN compilation flags

FFLAGS      = $(FOPT) -I$(ALICE_ROOT)/GEANT321

##### TARGETS #####
 
# Target

SLIBRARY	= $(ALICE_ROOT)/lib/libTGeant3.$(SL) \
		  $(ALICE_ROOT)/lib/libG3Dummy.$(SL) \
		  $(ALICE_ROOT)/lib/libTGeant3Dummy.$(SL) 
ALIBRARY	= $(ALICE_ROOT)/lib/libTGeant3.a

default:	headers $(SLIBRARY)

$(ALICE_ROOT)/lib/libTGeant3.$(SL):		$(GOBJS)

$(ALICE_ROOT)/lib/libG3Dummy.$(SL):		$(JOBJS)

$(ALICE_ROOT)/lib/libTGeant3Dummy.$(SL):	$(DOBJS)

$(GDICT):				$(GHDRS)

$(DDICT):				$(DHDRS)

DICT	=		$(GDICT) $(DDICT)

DEPINC	+= -I$(ALICE_ROOT)/GEANT321

HEADERS 	=	$(filter-out %LinkDef.h,$(GHDRS))

depend:			$(GSRCS) $(FSRCS)

TOCLEAN		= $(GOBJS) $(JOBJS) $(DOBJS) *Cint.cxx *Cint.h

############################### General Macros ################################

include $(ALICE_ROOT)/conf/GeneralMacros

############################ Dependencies #####################################

include Make-depend 
