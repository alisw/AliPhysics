############################### LHC Makefile ################################

# Include machine specific definitions

include $(ALICE_ROOT)/conf/GeneralDef
include $(ALICE_ROOT)/conf/MachineDef.$(ALICE_TARGET)

PACKAGE = LHC

# C++ sources

SRCS          = AliLHC.cxx AliLhcIRegion.cxx AliLhcBeam.cxx AliLhcProcess.cxx \
  	        AliLhcProcessBB.cxx AliLhcProcessIBS.cxx AliLhcMonitor.cxx \
	        AliLhcProcessBT.cxx

# C++ Headers

HDRS          = $(SRCS:.cxx=.h) 

# Library dictionary

DICT          = LHCCint.cxx
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

FFLAGS      = $(FOPT)

##### TARGETS #####
 
# Target

SLIBRARY	= $(LIBDIR)/libLHC.$(SL)
ALIBRARY	= $(LIBDIR)/libLHC.a

default:	$(SLIBRARY)

$(LIBDIR)/libLHC.$(SL):		$(OBJS)

$(DICT):				$(HDRS)

depend:					$(SRCS)


TOCLEAN		= $(OBJS) *Cint.cxx *Cint.h

CHECKS        = $(patsubst %.cxx,check/%.viol,$(SRCS))

############################### General Macros ################################

include $(ALICE_ROOT)/conf/GeneralMacros

############################ Dependencies #####################################

-include tgt_$(ALICE_TARGET)/Make-depend 







