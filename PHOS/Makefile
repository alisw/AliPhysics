############################### PHOS Makefile #################################

# Include machine specific definitions

include $(ALICE_ROOT)/conf/GeneralDef
include $(ALICE_ROOT)/conf/MachineDef.$(ALICE_TARGET)

PACKAGE = PHOS

# C++ sources

SRCS          =    AliPHOS.cxx AliPHOSv0.cxx AliPHOSv1.cxx AliPHOSv2.cxx

# C++ Headers

HDRS          = $(SRCS:.cxx=.h) PHOSLinkDef.h

# Library dictionary

DICT          = PHOSCint.cxx
DICTH         = $(DICT:.cxx=.h)
DICTO         = $(patsubst %.cxx,tgt_$(ALICE_TARGET)/%.o,$(DICT))

# FORTRAN sources

# FORTRAN Objects

# C Objects

# C++ Objects

OBJS          = $(patsubst %.cxx,tgt_$(ALICE_TARGET)/%.o,$(SRCS)) $(DICTO)

# C++ compilation flags

CXXFLAGS      = $(CXXOPTS) -I$(ROOTSYS)/include -I. -I$(ALICE_ROOT)/include/
 
# FORTRAN compilation flags

ALSRCS  = $(SRCS)

ALOBJS  = $(OBJS)

##### COMMANDS ##### 

SLIBRARY     	= $(LIBDIR)/libPHOS.$(SL)

default:	$(SLIBRARY)

$(LIBDIR)/libPHOS.$(SL):			$(OBJS) $(FOBJS)

$(DICT):					$(HDRS)

depend:				$(SRCS) 

TOCLEAN		=	$(OBJS) *Cint.cxx *Cint.h

############################### General Macros ################################

include $(ALICE_ROOT)/conf/GeneralMacros

############################ Dependencies #####################################

-include tgt_$(ALICE_TARGET)/Make-depend 

### Target check creates violation reports (.viol), which depend on
### stripped files (.ii), which in turn depend on preprocessed
### files (.i). Dependences are in conf/GeneralDef.

CHECKS        = $(patsubst %.cxx,check/%.viol,$(SRCS))

check:          $(CHECKS)







