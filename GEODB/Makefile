############################### GEODB Makefile ################################

# Include machine specific definitions

include $(ALICE_ROOT)/conf/GeneralDef
include $(ALICE_ROOT)/conf/MachineDef.$(ALICE_TARGET)

PACKAGE = GEODB

# C++ sources

SRCS          = AliGBox.cxx AliGCone.cxx AliGConfig.cxx \
         AliGMaterial.cxx AliGNode.cxx AliGShape.cxx AliGSphere.cxx \
         AliGSuperShape.cxx AliGTransform.cxx AliGTube.cxx AliGPCone.cxx \
         AliGTRD1.cxx AliGeometry.cxx AliGEODB.cxx


# C++ Headers

HDRS          = $(SRCS:.cxx=.h) GEODBLinkDef.h

# Library dictionary

DICT          = GEODBCint.cxx
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

SLIBRARY	= $(LIBDIR)/libGEODB.$(SL)
ALIBRARY	= $(LIBDIR)/libGEODB.a

default:	$(SLIBRARY) headers

$(LIBDIR)/libGEODB.$(SL):	$(OBJS)

$(DICT):			$(HDRS)

depend:				$(SRCS)

TOCLEAN			= $(OBJS) *Cint.cxx *Cint.h

CHECKS        = $(patsubst %.cxx,check/%.viol,$(SRCS))

############################### General Macros ################################

include $(ALICE_ROOT)/conf/GeneralMacros

############################ Dependencies #####################################

-include tgt_$(ALICE_TARGET)/Make-depend 
