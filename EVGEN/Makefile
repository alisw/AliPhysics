############################### EVGEN Makefile ################################

# Include machine specific definitions

include $(ALICE_ROOT)/conf/GeneralDef
include $(ALICE_ROOT)/conf/MachineDef.$(ALICE_TARGET)

PACKAGE = EVGEN

# C++ sources

SRCS          = AliSimpleGen.cxx AliGenCocktail.cxx AliGenFLUKAsource.cxx \
                AliGenMUONlib.cxx AliGenParam.cxx AliGenPythia.cxx \
                AliPythia.cxx AliDimuCombinator.cxx AliGenHalo.cxx \
                AliGenExtFile.cxx AliGenScan.cxx AliGenPHOSlib.cxx \
		AliGenDoubleScan.cxx
# C++ Headers

HDRS          = $(SRCS:.cxx=.h) EVGENLinkDef.h

# Library dictionary

DICT          = EVGENCint.cxx
DICTH         = $(DICT:.cxx=.h)
DICTO         = $(patsubst %.cxx,tgt_$(ALICE_TARGET)/%.o,$(DICT))

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

-include tgt_$(ALICE_TARGET)/Make-depend 

### Target check creates violation reports (.viol), which depend on
### stripped files (.ii), which in turn depend on preprocessed
### files (.i). Dependences are in conf/GeneralDef.

CHECKS        = $(patsubst %.cxx,check/%.viol,$(SRCS))

check:          $(CHECKS)







