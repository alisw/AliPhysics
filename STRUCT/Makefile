############################### STRUCT Makefile ###############################

# Include machine specific definitions

include $(ALICE_ROOT)/conf/GeneralDef
include $(ALICE_ROOT)/conf/MachineDef.$(ALICE_TARGET)

PACKAGE = STRUCT

# C++ sources

SRCS          = AliMAG.cxx      AliABSO.cxx       \
		AliDIPO.cxx  	AliDIPOv1.cxx     AliDIPOv2.cxx \
                AliFRAME.cxx    AliFRAMEv0.cxx    AliFRAMEv1.cxx \
                AliHALL.cxx       AliSHIL.cxx \
                AliPIPE.cxx     AliPIPEv0.cxx     AliPIPEv1.cxx \
		AliPIPEv3.cxx 	AliBODY.cxx

# C++ Headers

HDRS          = $(SRCS:.cxx=.h) STRUCTLinkDef.h

# Library dictionary

DICT          = STRUCTCint.cxx
DICTH         = $(DICT:.cxx=.h)
DICTO         = $(DICT:.cxx=.o)

# FORTRAN Objectrs

FOBJS         = $(FSRCS:.f=.o)

# C Objects

COBJS         = $(CSRCS:.c=.o)

# C++ Objects

OBJS          = $(SRCS:.cxx=.o) $(DICTO)

# C++ compilation flags

CXXFLAGS      = $(CXXOPTS) -I$(ROOTSYS)/include -I. -I$(ALICE_ROOT)/include/
 
# FORTRAN compilation flags

FFLAGS      = $(FOPT)

##### TARGETS #####
 
# Target

SLIBRARY	= $(ALICE_ROOT)/lib/libSTRUCT.$(SL)
ALIBRARY	= $(ALICE_ROOT)/lib/libSTRUCT.a

default:	$(SLIBRARY)

$(ALICE_ROOT)/lib/libSTRUCT.$(SL):		$(OBJS)

$(DICT):				$(HDRS)

depend:					$(SRCS)

TOCLEAN			= $(OBJS) *Cint.cxx *Cint.h

############################### General Macros ################################

include $(ALICE_ROOT)/conf/GeneralMacros

############################ Dependencies #####################################

-include Make-depend 
