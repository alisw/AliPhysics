############################### TPC Makefile ##################################

# Include machine specific definitions

include $(ALICE_ROOT)/conf/GeneralDef
include $(ALICE_ROOT)/conf/MachineDef.$(ALICE_TARGET)

PACKAGE = TPC

# C++ sources

SRCS          = AliTPC.cxx  AliTPCv0.cxx    AliTPCv1.cxx      AliTPCv2.cxx  \
		AliTPCv3.cxx AliDetectorParam.cxx AliTPCParam.cxx \
		AliTPCParamSR.cxx AliTPCParamCR.cxx  AliTPCPRF2D.cxx \
		AliTPCRF1D.cxx AliArrayI.cxx AliArrayS.cxx \
		AliSegmentID.cxx  AliSegmentArray.cxx AliDigits.cxx \
		AliSimDigits.cxx AliDigitsArray.cxx  AliTPCDigitsArray.cxx \
		AliCluster.cxx AliClusters.cxx AliClustersArray.cxx \
		AliClusterFinder.cxx AliTPCClustersRow.cxx \
		AliTPCClustersArray.cxx AliH2F.cxx


# C++ Headers

HDRS          = $(SRCS:.cxx=.h) TPCLinkDef.h

# Library dictionary

DICT          = TPCCint.cxx
DICTH         = $(DICT:.cxx=.h)
DICTO         = $(patsubst %.cxx,tgt_$(ALICE_TARGET)/%.o,$(DICT))

# FORTRAN Objectrs

FOBJS         = $(FSRCS:.f=.o)

# C Objects

COBJS         = $(CSRCS:.c=.o)

# C++ Objects

OBJS          = $(patsubst %.cxx,tgt_$(ALICE_TARGET)/%.o,$(SRCS)) $(DICTO)

# C++ compilation flags

CXXFLAGS      = $(CXXOPTS) -g -I$(ROOTSYS)/include -I. -I$(ALICE_ROOT)/include
#CXXFLAGS      = $(CXXOPTS) -I$(ROOTSYS)/include -I. -I$(ALICE_ROOT)/include/
 
# FORTRAN compilation flags

FFLAGS      = $(FOPT)
 
##### TARGETS #######
 
# Target

SLIBRARY	= $(LIBDIR)/libTPC.$(SL)

default:	$(SLIBRARY)

$(LIBDIR)/libTPC.$(SL):  $(OBJS)

$(DICT):		$(HDRS)

depend:			$(SRCS)

TOCLEAN			= $(OBJS) *Cint.h *Cint.cxx

############################### General Macros ################################

include $(ALICE_ROOT)/conf/GeneralMacros

############################ Dependencies #####################################

include tgt_$(ALICE_TARGET)/Make-depend 


