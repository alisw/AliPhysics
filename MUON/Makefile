############################### MUON Makefile #################################

# Include machine specific definitions

include $(ALICE_ROOT)/conf/GeneralDef
include $(ALICE_ROOT)/conf/MachineDef.$(ALICE_TARGET)

PACKAGE = MUON

# C++ sources

SRCS         = AliMUONChamber.cxx AliMUONChamberTrigger.cxx \
	       AliMUONSegmentationV0.cxx\
	       AliMUONResponse.cxx AliMUONResponseV0.cxx \
	       AliMUONSegmentationV01.cxx \
               AliMUONSegmentationV02.cxx AliMUONSegmentationV04.cxx \
	       AliMUONSegmentationV05.cxx\
	       AliMUONSegmentationTrigger.cxx  AliMUONResponseTrigger.cxx\
	       AliMUONSegmentationTriggerX.cxx AliMUONSegmentationTriggerY.cxx \
               AliMUONSegmentationV1.cxx AliMUON.cxx AliMUONv0.cxx AliMUONv1.cxx\
               AliMUONDisplay.cxx AliMUONPoints.cxx \
               AliMUONClusterFinderVS.cxx \
               AliMUONHitMapA1.cxx \
	       AliMUONHit.cxx AliMUONPadHit.cxx AliMUONDigit.cxx \
               AliMUONTransientDigit.cxx AliMUONRawCluster.cxx \
	       AliMUONReconstHit.cxx \
	       AliMUONTrackReconstructor.cxx \
               AliMUONEventReconstructor.cxx \
	       AliMUONTriggerDecision.cxx \
               AliMUONHitForRec.cxx AliMUONSegment.cxx \
               AliMUONTrack.cxx AliMUONTrackHit.cxx AliMUONTrackParam.cxx \
               AliMUONTriggerCircuit.cxx AliMUONTriggerLut.cxx \
               AliMUONGlobalTrigger.cxx AliMUONLocalTrigger.cxx \
               AliMUONTriggerConstants.cxx  AliMUONConstants.cxx \
	       AliMUONClusterInput.cxx     

# C++ Headers

HDRS          = $(SRCS:.cxx=.h) $(ROOTSYS)/include/TTree.h \
		$(ROOTSYS)/include/TMatrix.h MUONLinkDef.h

# Library dictionary

DICT          = MUONCint.cxx
DICTH         = $(DICT:.cxx=.h)
DICTO         = $(patsubst %.cxx,tgt_$(ALICE_TARGET)/%.o,$(DICT))

# Fortran sources

FSRCS	      = reco_muon.F extrap.F

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

SLIBRARY	= $(LIBDIR)/libMUON.$(SL)
ALIBRARY	= $(LIBDIR)/libMUON.a

default:	$(SLIBRARY)

$(LIBDIR)/libMUON.$(SL):	$(OBJS) $(FOBJS)

$(DICT):			$(HDRS)

depend:				$(SRCS)

TOCLEAN		= $(OBJS) $(FOBJS) *Cint.cxx *Cint.h

CHECKS        = $(patsubst %.cxx,check/%.viol,$(SRCS))

############################### General Macros ################################

include $(ALICE_ROOT)/conf/GeneralMacros

############################ Dependencies #####################################

-include tgt_$(ALICE_TARGET)/Make-depend 
