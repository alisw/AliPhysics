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
DICTO         = $(DICT:.cxx=.o)

# FORTRAN sources

FSRCS = AliPHOSf.F AliPHOScommon.F

# FORTRAN Objectrs

FOBJS         = $(FSRCS:.F=.o)

# C Objects

COBJS         = $(CSRCS:.c=.o)

# C++ Objects

OBJS          = $(SRCS:.cxx=.o) $(DICTO)

# C++ compilation flags

CXXFLAGS      = $(CXXOPTS) -I$(ROOTSYS)/include -I. -I$(ALICE_ROOT)/include/
 
# FORTRAN compilation flags

FFLAGS      = $(FOPT)

##### Other Libraries #####

SHSRCS  = $(wildcard shaker/*.f)
SHOBJS  = $(SHSRCS:.f=.o)

RCSRCS  = $(wildcard reconstruction/*.f) $(wildcard utils/*.cxx) $(wildcard utils/*.f) 
RCOBJS  = $(patsubst %.f,%.o,$(patsubst %.cxx,%.o,$(RCSRCS))) 

DUSRCS	= libPHOSdummy.cxx
DUOBJS	= $(DUSRCS:.cxx=.o)

ALSRCS  = $(SRCS) $(SHSRCS) $(RCSRCS) $(DUSRCS) dummies.c

ALOBJS  = $(SHOBJS) $(RCOBJS) $(DUOBJS)

##### COMMANDS ##### 

SLIBRARY     	= $(ALICE_ROOT)/lib/libPHOSdummy.$(SL) $(ALICE_ROOT)/lib/libPHOSshaker.$(SL) \
                  $(ALICE_ROOT)/lib/libPHOSreconstruction.$(SL) $(ALICE_ROOT)/lib/libPHOS.$(SL)

default:	$(SLIBRARY)

$(ALICE_ROOT)/lib/libPHOSdummy.$(SL):		$(DUOBJS)

$(ALICE_ROOT)/lib/libPHOSshaker.$(SL):		$(SHOBJS)

$(ALICE_ROOT)/lib/libPHOSreconstruction.$(SL):	$(RCOBJS)

$(ALICE_ROOT)/lib/libPHOS.$(SL):			$(OBJS) $(FOBJS)

$(DICT):					$(HDRS)

depend:				$(SRCS) $(DUSRCS) $(SHSRCS) $(RCSRCS)

TOCLEAN		=	$(DUOBJS) $(SHOBJS) $(RCOBJS) $(OBJS) $(FOBJS) \
                        *Cint.cxx *Cint.h

############################### General Macros ################################

include $(ALICE_ROOT)/conf/GeneralMacros

############################ Dependencies #####################################

include Make-depend 
