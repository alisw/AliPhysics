############################### EVGEN Makefile ################################

# Include machine specific definitions

include $(ALICE_ROOT)/conf/GeneralDef
include $(ALICE_ROOT)/conf/MachineDef.$(ALICE_TARGET)

PACKAGE = EVGEN

# C++ sources

SRCS          = AliSimpleGen.cxx AliGenCocktail.cxx AliGenFLUKAsource.cxx \
                AliGenMUONlib.cxx AliGenParam.cxx AliGenPythia.cxx \
                AliPythia.cxx AliDimuCombinator.cxx AliGenHalo.cxx \
                AliGenExtFile.cxx AliGenScan.cxx
# C++ Headers

HDRS          = $(SRCS:.cxx=.h) EVGENLinkDef.h

# Library dictionary

DICT          = EVGENCint.cxx
DICTH         = $(DICT:.cxx=.h)
DICTO         = $(patsubst %.cxx,$(ALICE_TARGET)/%.o,$(DICT))

# FORTRAN Objectrs

FOBJS         = $(patsubst %.F,$(ALICE_TARGET)/%.o,$(FSRCS))

# C Objects

COBJS         = $(patsubst %.c,$(ALICE_TARGET)/%.o,$(CSRCS))

# C++ Objects

OBJS          = $(patsubst %.cxx,$(ALICE_TARGET)/%.o,$(SRCS)) $(DICTO)

# C++ compilation flags

CXXFLAGS      = $(CXXOPTS) -I$(ROOTSYS)/include -I. -I$(ALICE_ROOT)/include/
 
# FORTRAN compilation flags

FFLAGS      = $(FOPT)

##### TARGETS #####
 
# Target

SLIBRARY	= $(LIBDIR)/libEVGEN.$(SL)
ALIBRARY	= $(LIBDIR)/libEVGEN.a

default:	$(SLIBRARY)

$(LIBDIR)/libEVGEN.$(SL):		$(OBJS)

$(DICT):				$(HDRS)

depend:					$(SRCS)

HEADERS		=	$(filter-out %LinkDef.h,$(HDRS))

TOCLEAN		= $(OBJS) *Cint.cxx *Cint.h

############################### General Macros ################################

include $(ALICE_ROOT)/conf/GeneralMacros

############################ Dependencies #####################################

-include Make-depend 
