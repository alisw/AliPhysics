############################### TOF Makefile ##################################

# Include machine specific definitions

include $(ALICE_ROOT)/conf/GeneralDef
include $(ALICE_ROOT)/conf/MachineDef.$(ALICE_TARGET)

PACKAGE = TOF

# C++ sources

SRCS  = AliTOF.cxx  AliTOFv0.cxx  AliTOFv1.cxx  AliTOFv2.cxx  AliTOFv3.cxx  AliTOFv4.cxx AliTOFv4T0.cxx AliTOFhit.cxx AliTOFhitT0.cxx AliTOFdigit.cxx AliTOFRawSector.cxx  AliTOFRoc.cxx  AliTOFRawDigit.cxx AliTOFDigitizer.cxx AliTOFSDigitizer.cxx AliTOFMerger.cxx  AliTOFv2FHoles.cxx AliTOFSDigit.cxx AliTOFHitMap.cxx AliTOFConstants.cxx AliTOFT0.cxx AliTOFPID.cxx AliTOFPad.cxx AliTOFRecHit.cxx AliTOFTrack.cxx AliTOFReconstructioner.cxx AliTOFProb.cxx

# C++ Headers

HDRS          = $(SRCS:.cxx=.h) TOFLinkDef.h

# Library dictionary

DICT          = TOFCint.cxx
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
 
# FORTRAN compilation flags

FFLAGS      = $(FOPT)

##### TARGETS #####
 
# Target

SLIBRARY	= $(LIBDIR)/libTOF.$(SL)
ALIBRARY	= $(LIBDIR)/libTOF.a

default:	$(SLIBRARY)

$(LIBDIR)/libTOF.$(SL):	$(OBJS)

$(DICT):			$(HDRS)

depend:				$(SRCS)

TOCLEAN			= $(OBJS) *Cint.cxx *Cint.h

CHECKS		        = $(patsubst %.cxx,check/%.viol,$(SRCS))

############################### General Macros ################################

include $(ALICE_ROOT)/conf/GeneralMacros

############################ Dependencies #####################################

-include tgt_$(ALICE_TARGET)/Make-depend 
# DO NOT DELETE
