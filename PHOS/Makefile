############################### PHOS Makefile #################################

# Include machine specific definitions

include $(ALICE_ROOT)/conf/GeneralDef
include $(ALICE_ROOT)/conf/MachineDef.$(ALICE_TARGET)

PACKAGE = PHOS

# C++ sources

SRCS          =    AliPHOS.cxx AliPHOSv0.cxx AliPHOSv1.cxx AliPHOSv2.cxx AliPHOSv3.cxx

# C++ Headers

HDRS          = $(SRCS:.cxx=.h) PHOSLinkDef.h

# Library dictionary

DICT          = PHOSCint.cxx
DICTH         = $(DICT:.cxx=.h)
DICTO         = $(patsubst %.cxx,tgt_$(ALICE_TARGET)/%.o,$(DICT))

# FORTRAN sources

FSRCS = AliPHOSf.F AliPHOScommon.F

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

##### Other Libraries #####

SHSRCS  = $(wildcard shaker/*.f)
SHOBJS  = $(patsubst %.f,tgt_$(ALICE_TARGET)/%.o,$(SHSRCS))

RCSRCS  = $(wildcard reconstruction/*.f)  
RCOBJS  = $(patsubst %.f,tgt_$(ALICE_TARGET)/%.o,$(RCSRCS))

DUSRCS	= libPHOSdummy.cxx
DUOBJS	= $(patsubst %.cxx,tgt_$(ALICE_TARGET)/%.o,$(DUSRCS))

ALSRCS  = $(SRCS) $(SHSRCS) $(RCSRCS) $(DUSRCS) dummies.c

ALOBJS  = $(SHOBJS) $(RCOBJS) $(DUOBJS)

##### COMMANDS ##### 

SLIBRARY     	= $(LIBDIR)/libPHOSdummy.$(SL) $(LIBDIR)/libPHOSshaker.$(SL) \
                  $(LIBDIR)/libPHOSreconstruction.$(SL) $(LIBDIR)/libPHOS.$(SL)

default:	$(SLIBRARY)

$(LIBDIR)/libPHOSdummy.$(SL):		$(DUOBJS)

$(LIBDIR)/libPHOSshaker.$(SL):		$(SHOBJS)

$(LIBDIR)/libPHOSreconstruction.$(SL):	$(RCOBJS)

$(LIBDIR)/libPHOS.$(SL):			$(OBJS) $(FOBJS)

$(DICT):					$(HDRS)

depend:				$(SRCS) $(DUSRCS) $(SHSRCS) $(RCSRCS)

TOCLEAN		=	$(DUOBJS) $(SHOBJS) $(RCOBJS) $(OBJS) $(FOBJS) \
                        *Cint.cxx *Cint.h

############################### General Macros ################################

include $(ALICE_ROOT)/conf/GeneralMacros

############################ Dependencies #####################################

-include Make-depend 
