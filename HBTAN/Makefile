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
                AliHBTQResolutionFctns.cxx AliHBTReaderPPprod.cxx\
		AliHBTReaderITSv1.cxx AliHBTReaderITSv2.cxx\
		AliHBTReaderKineTree.cxx AliHBTReaderInternal.cxx\
		AliHBTTwoTrackEffFctn.cxx \
		AliHBTLLWeights.cxx AliHBTLLWeightFctn.cxx

# Fortran sources

FSRCS	       = fsiini.F fsiw.F led_bldata.F ltran12.F

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

INCLUDES      = -I$(ALICE_ROOT)/TPC/ -I$(ALICE_ROOT)/CONTAINERS/ -I$(ALICE_ROOT)/ITS/

CXXFLAGS      = $(CXXOPTS) -I$(ROOTSYS)/include -I. -I$(ALICE_ROOT)/include/ $(INCLUDES)
 
# FORTRAN compilation flags

FFLAGS      = $(FOPT) -I$(ALICE_ROOT)/GEANT321

##### TARGETS #####
 
# Target

SLIBRARY	= $(LIBDIR)/libHBTAN.$(SL)
ALIBRARY	= $(LIBDIR)/libHBTAN.a

default:	$(SLIBRARY)


$(LIBDIR)/libHBTAN.$(SL):		$(OBJS)

$(DICT):			$(HDRS)

DEPINC += -I$(ALICE_ROOT)/GEANT321

depend:				$(SRCS) $(FSRCS)

TOCLEAN		= $(OBJS) $(FOBJS) *Cint.cxx *Cint.h

CHECKS        = $(patsubst %.cxx,check/%.viol,$(SRCS))

############################### General Macros ################################

include $(ALICE_ROOT)/conf/GeneralMacros

############################ Dependencies #####################################

-include tgt_$(ALICE_TARGET)/Make-depend 


