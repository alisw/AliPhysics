############################### ALIFAST Makefile ##################################

# Include machine specific definitions

include $(ALICE_ROOT)/conf/GeneralDef
include $(ALICE_ROOT)/conf/MachineDef.$(ALICE_TARGET)

PACKAGE = ALIFAST

# C++ sources

SRCS          = AliFast.cxx AliFMaker.cxx AliFDet.cxx AliFTrackMaker.cxx \
 		AliFTrack.cxx AliFHistBrowser.cxx AliFBigBang.cxx \
		AliFBrowsable.cxx AliFVirtualDisplay.cxx

# C++ Headers

HDRS          = $(SRCS:.cxx=.h) ALIFASTLinkDef.h

# Library dictionary

DICT          = ALIFASTCint.cxx
DICTH         = $(DICT:.cxx=.h)
DICTO         = $(patsubst %.cxx,tgt_$(ALICE_TARGET)/%.o,$(DICT))

# FORTRAN Objectrs

FOBJS         = $(FSRCS:.f=.o)

# C Objects

COBJS         = $(CSRCS:.c=.o)

# C++ Objects

OBJS          = $(patsubst %.cxx,tgt_$(ALICE_TARGET)/%.o,$(SRCS)) $(DICTO)

# C++ compilation flags

CXXFLAGS      = $(CXXOPTS) -I$(ROOTSYS)/include -I. -I$(ALICE_ROOT)/include/
 
# FORTRAN compilation flags

FFLAGS      = $(FOPT)
 
##### TARGETS #####
 
# Target

SLIBRARY	= $(LIBDIR)/libALIFAST.$(SL)

default:	$(SLIBRARY)

$(LIBDIR)/libALIFAST.$(SL):  $(OBJS)

$(DICT):		$(HDRS)

depend:			$(SRCS)

TOCLEAN			= $(OBJS) *Cint.h *Cint.cxx

############################### General Macros ################################

include $(ALICE_ROOT)/conf/GeneralMacros

CHECKS        = $(patsubst %.cxx,check/%.viol,$(SRCS))

############################ Dependencies #####################################

-include tgt_$(ALICE_TARGET)/Make-depend 


