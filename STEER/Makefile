############################### STEER Makefile ################################

# Include machine specific definitions

include $(ALICE_ROOT)/conf/GeneralDef
include $(ALICE_ROOT)/conf/MachineDef.$(ALICE_TARGET)

PACKAGE = STEER

# C++ sources

SRCS          = AliDetector.cxx  AliHeader.cxx    AliMagF.cxx GParticle.cxx \
                AliDigit.cxx     AliHit.cxx       AliPoints.cxx \
                AliDisplay.cxx   AliMC.cxx        AliRun.cxx AliGenerator.cxx \
		AliLego.cxx AliModule.cxx


# C++ Headers

HDRS          = $(SRCS:.cxx=.h) STEERLinkDef.h

# Library dictionary

DICT          = STEERCint.cxx
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

SLIBRARY	= $(ALICE_ROOT)/lib/libSTEER.$(SL)
ALIBRARY	= $(ALICE_ROOT)/lib/libSTEER.a

default:	headers $(SLIBRARY)

$(ALICE_ROOT)/lib/libSTEER.$(SL):		$(OBJS)

$(DICT):				$(HDRS)

depend:					$(SRCS)

HEADERS =	$(filter-out %LinkDef.h,$(HDRS) AliConst.h AliCallf77.h)

TOCLEAN		= $(OBJS) *Cint.cxx *Cint.h

############################### General Macros ################################

include $(ALICE_ROOT)/conf/GeneralMacros

############################ Dependencies #####################################

-include Make-depend 
