############################### ITS Makefile ##################################

# Include machine specific definitions

include $(ALICE_ROOT)/conf/GeneralDef
include $(ALICE_ROOT)/conf/MachineDef.$(ALICE_TARGET)

PACKAGE = ITS

# C++ sources

SRCS          = AliITS.cxx \
AliITSv1.cxx AliITSv3.cxx AliITSv5.cxx\
AliITShit.cxx \
AliITSgeom.cxx AliITSgeomSPD.cxx AliITSgeomSDD.cxx AliITSgeomSSD.cxx \
AliITSmodule.cxx AliITSsimulationFastPoints.cxx\
AliITSsimulation.cxx AliITSsimulationSPD.cxx \
AliITSsimulationSDD.cxx AliITSetfSDD.cxx\
AliITSsimulationSSD.cxx AliITSdcsSSD.cxx\
AliITSdictSSD.cxx \
AliITSdigit.cxx \
AliITSRawCluster.cxx \
AliITSstatistics.cxx AliITSstatistics2.cxx \
AliITSMap.cxx \
AliITSsegmentationSPD.cxx  AliITSsegmentationSDD.cxx AliITSsegmentationSSD.cxx\
AliITSresponseSPD.cxx  AliITSresponseSDD.cxx AliITSresponseSSD.cxx\
AliITSClusterFinder.cxx \
AliITSClusterFinderSPD.cxx AliITSClusterFinderSDD.cxx AliITSRawData.cxx \
AliITSHuffman.cxx \
AliITSClusterFinderSSD.cxx AliITSclusterSSD.cxx AliITSpackageSSD.cxx \
AliITSDetType.cxx \
AliITSgeomSPD300.cxx AliITSgeomSPD425.cxx
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

