############################### RICH Makefile #################################

# Include machine specific definitions

include $(ALICE_ROOT)/conf/GeneralDef
include $(ALICE_ROOT)/conf/MachineDef.$(ALICE_TARGET)

PACKAGE = RICH

# C++ sources

SRCS   = AliRICH.cxx AliRICHv0.cxx AliRICHv1.cxx AliRICHv2.cxx\
	 AliRICHDisplay.cxx AliRICHPoints.cxx AliRICHGeometry.cxx\
	 AliRICHSegmentationV0.cxx AliRICHTresholdMap.cxx\
	 AliRICHSegmentationV1.cxx AliRICHResponse.cxx AliRICHResponseV0.cxx \
	 AliRICHClusterFinder.cxx AliRICHHitMapA1.cxx\
	 AliRICHChamber.cxx AliRICHDetect.cxx AliRICHPatRec.cxx \
         AliRICHHit.cxx AliRICHCerenkov.cxx AliRICHPadHit.cxx \
	 AliRICHDigit.cxx AliRICHTransientDigit.cxx \
	 AliRICHRawCluster.cxx AliRICHRecHit1D.cxx AliRICHRecHit3D.cxx AliRICHEllipse.cxx 

# C++ Headers

#HDRS          = $(SRCS:.cxx=.h) $(ROOTSYS)/include/TClonesArray.h RICHLinkDef.h
HDRS          = $(SRCS:.cxx=.h)  RICHLinkDef.h

# Library dictionary

DICT          = RICHCint.cxx
DICTH         = $(DICT:.cxx=.h)
DICTO         = $(patsubst %.cxx,tgt_$(ALICE_TARGET)/%.o,$(DICT))

# FORTRAN Objectrs

FOBJS         = $(patsubst %.f,tgt_$(ALICE_TARGET)/%.o,$(FSRCS))

# C Objects

COBJS         = $(patsubst %.c,tgt_$(ALICE_TARGET)/%.o,$(CSRCS))

# C++ Objects

OBJS          = $(patsubst %.cxx,tgt_$(ALICE_TARGET)/%.o,$(SRCS)) $(DICTO)

# C++ compilation flags

CXXFLAGS      = $(CXXOPTS) -I$(ROOTSYS)/include -I. -I$(ALICE_ROOT)/include/

CXXFLAGS      += -I$(ALICE_ROOT)/TGeant3
DEPINC        += -I$(ALICE_ROOT)/TGeant3
 
# FORTRAN compilation flags

FFLAGS      = $(FOPT)

##### TARGETS #####
 
# Target

SLIBRARY	= $(LIBDIR)/libRICH.$(SL)
ALIBRARY	= $(LIBDIR)/libRICH.a

default:	$(SLIBRARY)

$(LIBDIR)/libRICH.$(SL):		$(OBJS)

$(DICT):				$(HDRS)

depend:					$(SRCS)

TOCLEAN		= $(OBJS) *Cint.cxx *Cint.h

CHECKS        = $(patsubst %.cxx,check/%.viol,$(SRCS))

############################### General Macros ################################

include $(ALICE_ROOT)/conf/GeneralMacros

############################ Dependencies #####################################

-include tgt_$(ALICE_TARGET)/Make-depend 






