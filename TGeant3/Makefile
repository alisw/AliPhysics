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
GDICTO         = $(patsubst %.cxx,$(ALICE_TARGET)/%.o,$(GDICT))

# Geant Objects

GOBJS          = $(patsubst %.cxx,$(ALICE_TARGET)/%.o,$(GSRCS)) $(FOBJS) $(GDICTO)

# FORTRAN Objectrs

FOBJS         = $(patsubst %.F,$(ALICE_TARGET)/%.o,$(FSRCS))

# C Objects

COBJS         = $(patsubst %.c,$(ALICE_TARGET)/%.o,$(CSRCS))

# C++ Objects

OBJS          = $(patsubst %.cxx,$(ALICE_TARGET)/%.o,$(SRCS)) $(DICTO)

# Dummy Geant321 class

# C++ sources

DSRCS          = TGeant3Dummy.cxx

# C++ Headers

DHDRS          = TGeant3.h TGeant3LinkDef.h

# Library dictionary

DDICT          = TGeant3DummyCint.cxx
DDICTH         = $(DDICT:.cxx=.h)
DDICTO         = $(patsubst %.cxx,$(ALICE_TARGET)/%.o,$(DDICT))

# Geant Objects

DOBJS          = $(patsubst %.cxx,$(ALICE_TARGET)/%.o,$(DSRCS)) $(DDICTO)

# Dummy Routines in Geant321

JSRCS		= G3Dummy.cxx
JOBJS		= $(patsubst %.cxx,$(ALICE_TARGET)/%.o,$(JSRCS))

# C++ compilation flags

CXXFLAGS      = $(CXXOPTS) -I$(ROOTSYS)/include -I. -I$(ALICE_ROOT)/include/
 
# FORTRAN compilation flags

FFLAGS      = $(FOPT) -I$(ALICE_ROOT)/GEANT321

##### TARGETS #####
 
# Target

SLIBRARY	= $(LIBDIR)/libTGeant3.$(SL) \
		  $(LIBDIR)/libG3Dummy.$(SL) \
		  $(LIBDIR)/libTGeant3Dummy.$(SL) 
ALIBRARY	= $(LIBDIR)/libTGeant3.a

default:	$(SLIBRARY)

$(LIBDIR)/libTGeant3.$(SL):		$(GOBJS)

$(LIBDIR)/libG3Dummy.$(SL):		$(JOBJS)

$(LIBDIR)/libTGeant3Dummy.$(SL):	$(DOBJS)

$(GDICT):				$(GHDRS)

$(DDICT):				$(DHDRS)

DICT	=		$(GDICT) $(DDICT)

DEPINC	+= -I$(ALICE_ROOT)/GEANT321

depend:			$(GSRCS) $(FSRCS)

TOCLEAN		= $(GOBJS) $(JOBJS) $(DOBJS) *Cint.cxx *Cint.h

############################### General Macros ################################

include $(ALICE_ROOT)/conf/GeneralMacros

############################ Dependencies #####################################

-include Make-depend 
