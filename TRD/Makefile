############################### TRD Makefile ##################################

# Include machine specific definitions

include $(ALICE_ROOT)/conf/GeneralDef
include $(ALICE_ROOT)/conf/MachineDef.$(ALICE_TARGET)

# Package name

PACKAGE	      = TRD

# C++ sources

SRCS          = AliTRD.cxx AliTRDv0.cxx AliTRDv1.cxx             \
                AliTRDhit.cxx AliTRDdigit.cxx                    \
		AliTRDpixel.cxx AliTRDmatrix.cxx                 \
                AliTRDgeometry.cxx AliTRDgeometryFull.cxx        \
                AliTRDgeometryHole.cxx AliTRDdigitizer.cxx       \
                AliTRDclusterizer.cxx AliTRDclusterizerV0.cxx    \
                AliTRDclusterizerV1.cxx AliTRDrecPoint.cxx       \
                AliTRDsegmentArray.cxx AliTRDdataArray.cxx       \
                AliTRDdataArrayI.cxx AliTRDdataArrayF.cxx        \
                AliTRDdigitsManager.cxx AliTRDsim.cxx            \
                AliTRDsegmentID.cxx AliTRDsegmentArrayBase.cxx   \
                AliTRDarrayI.cxx AliTRDarrayF.cxx                \
                AliTRDpoints.cxx AliTRDtimeBin.cxx               \
                AliTRDtrackingSector.cxx                         \
                AliTRDtrack.cxx AliTRDtracker.cxx                \
                AliTRDcluster.cxx         

# C++ Headers

HDRS          = $(SRCS:.cxx=.h)                \
                $(ROOTSYS)/include/TF1.h       \
                $(ROOTSYS)/include/TFile.h     \
                $(ROOTSYS)/include/TObjArray.h \
                TRDLinkDef.h

# Library dictionary

DICT          = TRDCint.cxx
DICTH         = $(DICT:.cxx=.h)
DICTO         = $(patsubst %.cxx,tgt_$(ALICE_TARGET)/%.o,$(DICT))

# FORTRAN Objectrs

FOBJS         = $(patsubst %.f,tgt_$(ALICE_TARGET)/%.o,$(FSRCS))

# C Objects

COBJS         = $(patsubst %.c,tgt_$(ALICE_TARGET)/%.o,$(CSRCS))

# C++ Objects

OBJS          = $(patsubst %.cxx,tgt_$(ALICE_TARGET)/%.o,$(SRCS)) $(DICTO)

# C++ compilation flags

CXXFLAGS      = $(CXXOPTS) -I$(ROOTSYS)/include -I. -I$(ALICE_ROOT)/include
 
# FORTRAN compilation flags

FFLAGS      = $(FOPT)

##### TARGETS #####
 
# Target

SLIBRARY	= $(LIBDIR)/libTRD.$(SL)
ALIBRARY	= $(LIBDIR)/libTRD.a

default:	$(SLIBRARY)

$(LIBDIR)/libTRD.$(SL):	$(OBJS)

TRDCint.cxx:			$(HDRS)

depend:				$(SRCS)

TOCLEAN			= $(OBJS) *Cint.h *Cint.cxx

CHECKS        = $(patsubst %.cxx,check/%.viol,$(SRCS))

############################### General Macros ################################

include $(ALICE_ROOT)/conf/GeneralMacros

############################ Dependencies #####################################

-include tgt_$(ALICE_TARGET)/Make-depend 








