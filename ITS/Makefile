############################### ITS Makefile ##################################

# Include machine specific definitions

include $(ALICE_ROOT)/conf/GeneralDef
include $(ALICE_ROOT)/conf/MachineDef.$(ALICE_TARGET)

PACKAGE = ITS

# C++ sources

SRCS          = AliITS.cxx AliITSv0.cxx    AliITSv1.cxx      AliITSv3.cxx    AliITSv4.cxx AliITSv5.cxx  

# Fortran sources

FSRCS	       = AliITSf.F 

# C++ Headers

HDRS          = $(SRCS:.cxx=.h) ITSLinkDef.h

# Library dictionary

DICT          = ITSCint.cxx
DICTH         = $(DICT:.cxx=.h)
DICTO         = $(DICT:.cxx=.o)

# FORTRAN Objectrs

FOBJS         = $(FSRCS:.F=.o)

# C Objects

COBJS         = $(CSRCS:.c=.o)

# C++ Objects

OBJS          = $(SRCS:.cxx=.o) $(DICTO)

# C++ compilation flags

CXXFLAGS      = $(CXXOPTS) -I$(ROOTSYS)/include -I. -I$(ALICE_ROOT)/include/
 
# FORTRAN compilation flags

FFLAGS      = $(FOPT) -I$(ALICE_ROOT)/GEANT321

##### TARGETS #####
 
# Target

SLIBRARY	= $(ALICE_ROOT)/lib/libITS.$(SL)
ALIBRARY	= $(ALICE_ROOT)/lib/libITS.a

default:	$(SLIBRARY)

$(ALICE_ROOT)/lib/libITS.$(SL):	$(OBJS) $(FOBJS)

$(DICT):			$(HDRS)

DEPINC += -I$(ALICE_ROOT)/GEANT321

depend:				$(SRCS) $(FSRCS)

TOCLEAN		= $(OBJS) $(FOBJS) *Cint.cxx *Cint.h

############################### General Macros ################################

include $(ALICE_ROOT)/conf/GeneralMacros

############################ Dependencies #####################################

include Make-depend 
