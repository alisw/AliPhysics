############################### STEER Makefile ################################

# Include machine specific definitions

include $(ALICE_ROOT)/conf/GeneralDef
include $(ALICE_ROOT)/conf/MachineDef.$(ALICE_TARGET)

PACKAGE = STEER

# C++ sources

SRCS          = AliDetector.cxx	AliHeader.cxx	AliMagF.cxx \
                AliDigit.cxx	AliHit.cxx	AliPoints.cxx \
                AliDisplay.cxx	AliMC.cxx	AliRun.cxx AliGenerator.cxx \
		AliLego.cxx	AliModule.cxx	AliDigitNew.cxx \
		AliGeometry.cxx	AliRecPoint.cxx	AliSegmentation.cxx \
	        AliHitMap.cxx	AliMagFC.cxx	AliMagFCM.cxx \
		AliMagFDM.cxx	AliLegoGenerator.cxx AliLegoGeneratorXYZ.cxx\
                AliLegoGeneratorEta.cxx 	AliRndm.cxx \
		AliKalmanTrack.cxx		AliCluster.cxx \
		AliMCQA.cxx

# C++ Headers

HDRS          = $(SRCS:.cxx=.h) AliPDG.h $(ROOTSYS)/include/TH2.h \
				$(ROOTSYS)/include/TRandom.h STEERLinkDef.h 

# Library dictionary

DICT          = STEERCint.cxx
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

SLIBRARY	= $(LIBDIR)/libSTEER.$(SL)
ALIBRARY	= $(LIBDIR)/lib/libSTEER.a

default:	$(SLIBRARY)

$(LIBDIR)/libSTEER.$(SL):		$(OBJS)

$(DICT):				$(HDRS)

depend:					$(SRCS)

HEADERS 	= $(SRCS:.cxx=.h) AliPDG.h AliConst.h AliCallf77.h AliMCProcess.h

TOCLEAN		= $(OBJS) *Cint.cxx *Cint.h

CHECKS        	= $(patsubst %.cxx,check/%.viol,$(SRCS))

############################### General Macros ################################

include $(ALICE_ROOT)/conf/GeneralMacros

############################ Dependencies #####################################

-include tgt_$(ALICE_TARGET)/Make-depend 







