############################### EVGEN Makefile ################################

# Include machine specific definitions

include $(ALICE_ROOT)/conf/GeneralDef
include $(ALICE_ROOT)/conf/MachineDef.$(ALICE_TARGET)

PACKAGE = EVGEN

# C++ sources

SRCS          = AliSimpleGen.cxx AliGenCocktail.cxx AliGenFLUKAsource.cxx \
                AliGenMUONlib.cxx AliGenParam.cxx AliGenPythia.cxx \
                AliPythia.cxx DimuonCombinator.cxx AliGenHalo.cxx \
                AliGenExtFile.cxx AliGenScan.cxx
# C++ Headers

HDRS          = $(SRCS:.cxx=.h) EVGENLinkDef.h

# Library dictionary

DICT          = EVGENCint.cxx
DICTH         = $(DICT:.cxx=.h)
DICTO         = $(DICT:.cxx=.o)

# FORTRAN Objectrs

FOBJS         = $(FSRCS:.F=.o)

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

SLIBRARY	= $(ALICE_ROOT)/lib/libEVGEN.$(SL)
ALIBRARY	= $(ALICE_ROOT)/lib/libEVGEN.a

default:	$(SLIBRARY) headers

$(ALICE_ROOT)/lib/libEVGEN.$(SL):		$(OBJS)

$(DICT):				$(HDRS)

depend:					$(SRCS)

HEADERS		=	$(filter-out %LinkDef.h,$(HDRS))

TOCLEAN		= $(OBJS) *Cint.cxx *Cint.h

############################### General Macros ################################

include $(ALICE_ROOT)/conf/GeneralMacros

############################ Dependencies #####################################

include Make-depend 
