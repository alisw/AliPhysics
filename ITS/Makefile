############################### ITS Makefile ##################################

# Include machine specific definitions

include $(ALICE_ROOT)/conf/GeneralDef
include $(ALICE_ROOT)/conf/MachineDef.$(ALICE_TARGET)

PACKAGE = ITS

# C++ sources

SRCS          = AliITS.cxx AliITSv1.cxx AliITSv3.cxx AliITSv5.cxx \
      AliITSv5symm.cxx AliITSv5asymm.cxx \
      AliITSvPPRcoarseasymm.cxx AliITSvPPRcoarsesymm.cxx\
      AliITShit.cxx AliITSDetType.cxx \
		AliITSgeom.cxx AliITSgeomMatrix.cxx \
      AliITSgeomSPD.cxx AliITSgeomSDD.cxx \
		AliITSgeomSSD.cxx AliITSmodule.cxx \
		AliITSsimulationFastPoints.cxx \
		AliITSsimulationFastPointsV0.cxx AliITSsimulation.cxx \
		AliITSsimulationSPD.cxx AliITSsimulationSPDbari.cxx \
		AliITSsimulationSDD.cxx \
		AliITSetfSDD.cxx AliITSsimulationSSD.cxx AliITSdcsSSD.cxx \
		AliITSdigit.cxx AliITSRawCluster.cxx AliITSRecPoint.cxx \
		AliITSMap.cxx AliITSMapA1.cxx AliITSMapA2.cxx \
		AliITSsegmentation.cxx AliITSsegmentationSPD.cxx \
		AliITSsegmentationSDD.cxx AliITSsegmentationSSD.cxx\
		AliITSresponse.cxx AliITSresponseSPD.cxx \
		AliITSresponseSPDbari.cxx \
		AliITSresponseSDD.cxx AliITSresponseSSD.cxx \
		AliITSClusterFinder.cxx AliITSClusterFinderSPD.cxx \
		AliITSClusterFinderSPDbari.cxx \
		AliITSClusterFinderSDD.cxx AliITSRawData.cxx \
		AliITSHuffman.cxx AliITSClusterFinderSSD.cxx \
		AliITSclusterSSD.cxx AliITSpackageSSD.cxx \
		AliITSdictSSD.cxx AliITSgeomSPD300.cxx AliITSgeomSPD425.cxx \
		AliITSstatistics.cxx AliITSstatistics2.cxx \
                AliITStrack.cxx AliITStracking.cxx AliITSiotrack.cxx		
#		AliITSAlignmentTrack.cxx AliITSAlignmentModule.cxx \
#		AliITSvtest.cxx AliITStrack.cxx

# Fortran sources

FSRCS	       = 

# C++ Headers

HDRS          = $(SRCS:.cxx=.h) ITSLinkDef.h

# Library dictionary

DICT          = ITSCint.cxx
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

FFLAGS      = $(FOPT) -I$(ALICE_ROOT)/GEANT321

##### TARGETS #####
 
# Target

SLIBRARY	= $(LIBDIR)/libITS.$(SL)
ALIBRARY	= $(LIBDIR)/libITS.a

default:	$(SLIBRARY)


$(LIBDIR)/libITS.$(SL):		$(OBJS)

$(DICT):			$(HDRS)

DEPINC += -I$(ALICE_ROOT)/GEANT321

depend:				$(SRCS) $(FSRCS)

TOCLEAN		= $(OBJS) $(FOBJS) *Cint.cxx *Cint.h

CHECKS        = $(patsubst %.cxx,check/%.viol,$(SRCS))

############################### General Macros ################################

include $(ALICE_ROOT)/conf/GeneralMacros

############################ Dependencies #####################################

-include tgt_$(ALICE_TARGET)/Make-depend 


