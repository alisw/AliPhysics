############################### HBTAnalysis Makefile ##################################

# Include machine specific definitions

include $(ALICE_ROOT)/conf/GeneralDef
include $(ALICE_ROOT)/conf/MachineDef.$(ALICE_TARGET)


OPT = -g
PACKAGE = HBTAnalysis

# C++ sources

SRCS          = AliHBTAnalysis.cxx AliHBTFunction.cxx \
                AliHBTEvent.cxx AliHBTRun.cxx \
		AliHBTParticle.cxx AliHBTParticleCut.cxx \
		AliHBTPair.cxx AliHBTPairCut.cxx\
                AliHBTCorrelFctn.cxx \
                AliHBTReader.cxx AliHBTReaderTPC.cxx\
                AliHBTQResolutionFctns.cxx AliHBTReaderPPprod.cxx

# Fortran sources

FSRCS	       = 

# C++ Headers

HDRS          = $(SRCS:.cxx=.h) HBTAnalysisLinkDef.h

# Library dictionary

DICT          = HBTAnalysisCint.cxx
DICTH         = $(DICT:.cxx=.h)
DICTO         = $(patsubst %.cxx,tgt_$(ALICE_TARGET)/%.o,$(DICT))

# FORTRAN Objectrs

FOBJS         = $(patsubst %.F,tgt_$(ALICE_TARGET)/%.o,$(FSRCS))

# C Objects

COBJS         = $(patsubst %.c,tgt_$(ALICE_TARGET)/%.o,$(CSRCS))

# C++ Objects

OBJS          = $(patsubst %.cxx,tgt_$(ALICE_TARGET)/%.o,$(SRCS)) $(DICTO)

# C++ compilation flags

INCLUDES      = -I$(ALICE_ROOT)/TPC/ -I$(ALICE_ROOT)/CONTAINERS/

CXXFLAGS      = $(CXXOPTS) -I$(ROOTSYS)/include -I. -I$(ALICE_ROOT)/include/ $(INCLUDES)
 
# FORTRAN compilation flags

FFLAGS      = $(FOPT) -I$(ALICE_ROOT)/GEANT321

##### TARGETS #####
 
# Target

SLIBRARY	= $(LIBDIR)/libHBTAnalysis.$(SL)
ALIBRARY	= $(LIBDIR)/libHBTAnalysis.a

default:	$(SLIBRARY)


$(LIBDIR)/libHBTAnalysis.$(SL):		$(OBJS)

$(DICT):			$(HDRS)

DEPINC += -I$(ALICE_ROOT)/GEANT321

depend:				$(SRCS) $(FSRCS)

TOCLEAN		= $(OBJS) $(FOBJS) *Cint.cxx *Cint.h

CHECKS        = $(patsubst %.cxx,check/%.viol,$(SRCS))

############################### General Macros ################################

include $(ALICE_ROOT)/conf/GeneralMacros

############################ Dependencies #####################################

-include tgt_$(ALICE_TARGET)/Make-depend 


