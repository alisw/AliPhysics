############################### EVGEN Makefile ################################

# Include machine specific definitions

include $(ALICE_ROOT)/conf/GeneralDef
include $(ALICE_ROOT)/conf/MachineDef.$(ALICE_TARGET)

PACKAGE = EVGEN

# C++ sources

SRCS          = AliGenHIJINGpara.cxx AliGenBox.cxx AliGenFixed.cxx \
		AliGenCocktail.cxx AliGenFLUKAsource.cxx \
		AliGenLib.cxx \
                AliGenMUONlib.cxx AliGenParam.cxx AliGenPythia.cxx \
                AliPythia.cxx AliDimuCombinator.cxx AliGenHalo.cxx \
                AliGenExtFile.cxx AliGenScan.cxx AliGenPHOSlib.cxx \
		AliGenDoubleScan.cxx AliGenCocktailEntry.cxx \
		AliGenGSIlib.cxx AliGenPMDlib.cxx\
	        AliGenHijing.cxx \
		AliGenEventHeader.cxx AliGenHijingEventHeader.cxx \
		AliDecayer.cxx AliDecayerPythia.cxx
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

HEADERS       = AliPythia.h GenTypeDefs.h AliDecayer.h AliDecayerPythia.h

TOCLEAN		= $(OBJS) *Cint.cxx *Cint.h

CHECKS        = $(patsubst %.cxx,check/%.viol,$(SRCS))

############################### General Macros ################################

include $(ALICE_ROOT)/conf/GeneralMacros

############################ Dependencies #####################################

-include tgt_$(ALICE_TARGET)/Make-depend 







