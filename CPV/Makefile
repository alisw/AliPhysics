############################### PHOS Makefile #################################

# Include machine specific definitions

include $(ALICE_ROOT)/conf/GeneralDef
include $(ALICE_ROOT)/conf/MachineDef.$(ALICE_TARGET)

PACKAGE = CPV

# C++ sources

SRCS          =    AliCPV.cxx AliCPVv0.cxx

# C++ Headers

HDRS          = $(SRCS:.cxx=.h) CPVLinkDef.h

# Library dictionary

DICT          = CPVCint.cxx
DICTH         = $(DICT:.cxx=.h)
DICTO         = $(patsubst %.cxx,tgt_$(ALICE_TARGET)/%.o,$(DICT))

# FORTRAN sources

FSRCS = cpvrec.F

# FORTRAN Objects

FOBJS         = $(patsubst %.F,tgt_$(ALICE_TARGET)/%.o,$(FSRCS))

# C Objects

COBJS         = $(patsubst %.c,tgt_$(ALICE_TARGET)/%.o,$(CSRCS))

# C++ Objects

OBJS          = $(patsubst %.cxx,tgt_$(ALICE_TARGET)/%.o,$(SRCS)) $(DICTO)

# C++ compilation flags

CXXFLAGS      = $(CXXOPTS) -I$(ROOTSYS)/include -I. -I$(ALICE_ROOT)/include/ \
                           -I$(ALICE_ROOT)/RALICE/
 
# FORTRAN compilation flags

FFLAGS        = $(FOPT)

##### COMMANDS ##### 

SLIBRARY     	= $(LIBDIR)/libCPV.$(SL)

default:	$(SLIBRARY)

$(LIBDIR)/libCPV.$(SL):		$(OBJS) $(FOBJS)

$(DICT):			$(HDRS)

depend:				$(SRCS)

TOCLEAN		=	$(OBJS) $(FOBJS) \
                        *Cint.cxx *Cint.h

############################### General Macros ################################

include $(ALICE_ROOT)/conf/GeneralMacros

############################ Dependencies #####################################

-include tgt_$(ALICE_TARGET)/Make-depend 

### Target check creates violation reports (.viol), which depend on
### stripped files (.ii), which in turn depend on preprocessed
### files (.i). Dependences are in conf/GeneralDef.

CHECKS        = $(patsubst %.cxx,check/%.viol,$(SRCS))

check:          $(CHECKS)

