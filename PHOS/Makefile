############################### PHOS Makefile #################################

# Include machine specific definitions

include $(ALICE_ROOT)/conf/GeneralDef
include $(ALICE_ROOT)/conf/MachineDef.$(ALICE_TARGET)

PACKAGE = PHOS

# C++ sources

SRCS          =  AliPHOS.cxx AliPHOSv0.cxx AliPHOSv1.cxx AliPHOSv2.cxx \
                 AliPHOSv3.cxx AliPHOSv4.cxx \
                 AliPHOSHit.cxx \
                 AliPHOSGeometry.cxx \
		 AliPHOSEMCAGeometry.cxx \
		 AliPHOSCPVGeometry.cxx AliPHOSCPVBaseGeometry.cxx \
		 AliPHOSPPSDGeometry.cxx AliPHOSSupportGeometry.cxx \
		 AliPHOSCPVDigit.cxx AliPHOSDigit.cxx  \
                 AliPHOSRecPoint.cxx AliPHOSEmcRecPoint.cxx \
		 AliPHOSPpsdRecPoint.cxx AliPHOSCpvRecPoint.cxx \
		 AliPHOSClusterizer.cxx AliPHOSClusterizerv1.cxx  AliPHOSLink.cxx \
		 AliPHOSReconstructioner.cxx  AliPHOSTrackSegment.cxx \
                 AliPHOSTrackSegmentMaker.cxx AliPHOSTrackSegmentMakerv1.cxx \
                 AliPHOSRecParticle.cxx AliPHOSPID.cxx AliPHOSPIDv1.cxx \
                 AliPHOSAnalyze.cxx AliPHOSFastRecParticle.cxx \
                 AliPHOSIndexToObject.cxx

# C++ Headers

HDRS          = $(SRCS:.cxx=.h) $(ROOTSYS)/include/TTree.h PHOSLinkDef.h

# Library dictionary

DICT          = PHOSCint.cxx
DICTH         = $(DICT:.cxx=.h)
DICTO         = $(patsubst %.cxx,tgt_$(ALICE_TARGET)/%.o,$(DICT))

# C++ Objects

OBJS          = $(patsubst %.cxx,tgt_$(ALICE_TARGET)/%.o,$(SRCS)) $(DICTO)

# C++ compilation flags

CXXFLAGS      = $(CXXOPTS) -I$(ROOTSYS)/include -I. -I$(ALICE_ROOT)/include/

ALSRCS  = $(SRCS) $(SHSRCS) $(RCSRCS) $(DUSRCS) dummies.c

ALOBJS  = $(SHOBJS) $(RCOBJS) $(DUOBJS)

##### COMMANDS ##### 

SLIBRARY     	= $(LIBDIR)/libPHOS.$(SL)

default:	$(SLIBRARY)

$(LIBDIR)/libPHOS.$(SL):			$(OBJS) $(FOBJS)

$(DICT):					$(HDRS)

depend:				$(SRCS) 

TOCLEAN		=	 $(OBJS) *Cint.cxx *Cint.h

CHECKS        = $(patsubst %.cxx,check/%.viol,$(SRCS))
REVENGS       = $(patsubst %.cxx,check/%.ii,$(SRCS))

# IRST reverse engineering
reveng : $(REVENGS)
	@$(REV_ENG) $(REVENGS)
############################### General Macros ################################

include $(ALICE_ROOT)/conf/GeneralMacros

############################ Dependencies #####################################

-include tgt_$(ALICE_TARGET)/Make-depend 







