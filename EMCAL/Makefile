############################### EMCAL Makefile #################################

# Include machine specific definitions

include $(ALICE_ROOT)/conf/GeneralDef
include $(ALICE_ROOT)/conf/MachineDef.$(ALICE_TARGET)

PACKAGE = EMCAL

# C++ sources

SRCS          =  AliEMCAL.cxx  AliEMCALGeometry.cxx AliEMCALv0.cxx 

# C++ Headers

HDRS          = $(SRCS:.cxx=.h) $(ROOTSYS)/include/TTree.h EMCALLinkDef.h

# Library dictionary

DICT          = EMCALCint.cxx
DICTH         = $(DICT:.cxx=.h)
DICTO         = $(patsubst %.cxx,tgt_$(ALICE_TARGET)/%.o,$(DICT))

# C++ Objects

OBJS          = $(patsubst %.cxx,tgt_$(ALICE_TARGET)/%.o,$(SRCS)) $(DICTO)

# C++ compilation flags

CXXFLAGS      = $(CXXOPTS) -I$(ROOTSYS)/include -I. -I$(ALICE_ROOT)/include/

ALSRCS  = $(SRCS) $(SHSRCS) $(RCSRCS) $(DUSRCS) dummies.c

ALOBJS  = $(SHOBJS) $(RCOBJS) $(DUOBJS)

##### COMMANDS ##### 

SLIBRARY     	= $(LIBDIR)/libEMCAL.$(SL)

default:	$(SLIBRARY)

$(LIBDIR)/libEMCAL.$(SL):			$(OBJS) $(FOBJS)

$(DICT):					$(HDRS)

depend:				$(SRCS) 

TOCLEAN		=	 $(OBJS) *Cint.cxx *Cint.h

CHECKS        = $(patsubst %.cxx,check/%.viol,$(SRCS))

############################### General Macros ################################

include $(ALICE_ROOT)/conf/GeneralMacros

############################ Dependencies #####################################

-include tgt_$(ALICE_TARGET)/Make-depend 







