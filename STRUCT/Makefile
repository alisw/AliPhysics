############################### STRUCT Makefile ###############################

# Include machine specific definitions

include $(ALICE_ROOT)/conf/GeneralDef
include $(ALICE_ROOT)/conf/MachineDef.$(ALICE_TARGET)

PACKAGE = STRUCT

# C++ sources

SRCS          = AliMAG.cxx      AliABSO.cxx       AliABSOv0.cxx \
		AliDIPO.cxx  	AliDIPOv1.cxx     AliDIPOv2.cxx \
                AliFRAME.cxx    AliFRAMEv0.cxx    AliFRAMEv1.cxx \
                AliHALL.cxx     AliSHIL.cxx 	  AliSHILv0.cxx AliSHILvF.cxx\
                AliPIPE.cxx     AliPIPEv0.cxx     AliPIPEv1.cxx \
		AliPIPEv3.cxx 	AliBODY.cxx	  AliALIFE.cxx

# C++ Headers

HDRS          = $(SRCS:.cxx=.h) STRUCTLinkDef.h

# Library dictionary

DICT          = STRUCTCint.cxx
DICTH         = $(DICT:.cxx=.h)
DICTO         = $(patsubst %.cxx,tgt_$(ALICE_TARGET)/%.o,$(DICT))

# FORTRAN Objectrs

FOBJS         = $(patsubst %.f,tgt_$(ALICE_TARGET)/%.o,$(FSRCS))

# C Objects

COBJS         = $(patsubst %.cxx,tgt_$(ALICE_TARGET)/%.o,$(CSRCS))

# C++ Objects

OBJS          = $(patsubst %.cxx,tgt_$(ALICE_TARGET)/%.o,$(SRCS)) $(DICTO)

# C++ compilation flags

CXXFLAGS      = $(CXXOPTS) -I$(ROOTSYS)/include -I. -I$(ALICE_ROOT)/include/
 
# FORTRAN compilation flags

FFLAGS      = $(FOPT)

##### TARGETS #####
 
# Target

SLIBRARY	= $(LIBDIR)/libSTRUCT.$(SL)
ALIBRARY	= $(LIBDIR)/libSTRUCT.a

default:	$(SLIBRARY)

$(LIBDIR)/libSTRUCT.$(SL):		$(OBJS)

$(DICT):				$(HDRS)

depend:					$(SRCS)

TOCLEAN			= $(OBJS) *Cint.cxx *Cint.h

CHECKS        = $(patsubst %.cxx,check/%.viol,$(SRCS))

############################### General Macros ################################

include $(ALICE_ROOT)/conf/GeneralMacros

############################ Dependencies #####################################

-include tgt_$(ALICE_TARGET)/Make-depend 







