############################### TRD Makefile ##################################

# Include machine specific definitions

include $(ALICE_ROOT)/conf/GeneralDef
include $(ALICE_ROOT)/conf/MachineDef.$(ALICE_TARGET)

# Package name

PACKAGE	      = TRD

# C++ sources

SRCS          = AliTRD.cxx AliTRDv0.cxx AliTRDv1.cxx           \
		AliTRDpixel.cxx AliTRDmatrix.cxx               \
                AliTRDgeometry.cxx AliTRDgeometryFull.cxx      \
                AliTRDgeometryHole.cxx AliTRDdigitizer.cxx     \
                AliTRDclusterizer.cxx AliTRDclusterizerV0.cxx  \
                AliTRDclusterizerV1.cxx AliTRDrecPoint.cxx     \
                AliTRDsegmentArray.cxx AliTRDdataArray.cxx     \
                AliTRDsegmentID.cxx AliTRDsegmentArrayBase.cxx \
                AliTRDarrayI.cxx

# C++ Headers

HDRS          = $(SRCS:.cxx=.h) TRDLinkDef.h

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

############################### General Macros ################################

include $(ALICE_ROOT)/conf/GeneralMacros

############################ Dependencies #####################################

-include tgt_$(ALICE_TARGET)/Make-depend 

### Target check creates violation reports (.viol), which depend on
### stripped files (.ii), which in turn depend on preprocessed
### files (.i). Dependences are in conf/GeneralDef.

CHECKS        = $(patsubst %.cxx,check/%.viol,$(SRCS))

check:          $(CHECKS)







