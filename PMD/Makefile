############################### PMD Makefile ##################################

# Include machine specific definitions

include $(ALICE_ROOT)/conf/GeneralDef
include $(ALICE_ROOT)/conf/MachineDef.$(ALICE_TARGET)

PACKAGE = PMD

# C++ sources

SRCS          = AliPMD.cxx AliPMDv0.cxx AliPMDv1.cxx AliPMDv2.cxx AliPMDRecPoint.cxx

# C++ Headers

HDRS    = $(SRCS:.cxx=.h) \
	$(ROOTSYS)/include/TClonesArray.h \
	PMDLinkDef.h

# Library dictionary

DICT          = PMDCint.cxx
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

SLIBRARY	= $(LIBDIR)/libPMD.$(SL)
ALIBRARY	= $(LIBDIR)/libPMD.a

default:	$(SLIBRARY)

$(LIBDIR)/libPMD.$(SL):	$(OBJS)

$(DICT):			$(HDRS)

depend:				$(SRCS)

TOCLEAN		= $(OBJS) *Cint.cxx *Cint.h

CHECKS        = $(patsubst %.cxx,check/%.viol,$(SRCS))

############################### General Macros ################################

include $(ALICE_ROOT)/conf/GeneralMacros

############################ Dependencies #####################################

-include tgt_$(ALICE_TARGET)/Make-depend 







