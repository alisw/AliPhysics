############################### AliGeant3 Makefile ##############################

# Include machine specific definitions

include $(ALICE_ROOT)/conf/GeneralDef
include $(ALICE_ROOT)/conf/MachineDef.$(ALICE_TARGET)

PACKAGE = AliGeant3

# FORTRAN sources

FSRCS	     = 

# C++ sources

GSRCS          = AliG3Medium.cxx \
                 AliG3Material.cxx AliG3Volume.cxx \
		 AliGUISliders.cxx AliGuiGeomDialog.cxx AliGuiGeomMain.cxx \
		 AliGeant3GeometryGUI.cxx AliNode.cxx AliG3toRoot.cxx

# C++ Headers

GHDRS          = $(filter-out gucode.h,$(GSRCS:.cxx=.h)) AliGeant3LinkDef.h

# Library dictionary

GDICT          = AliGeant3Cint.cxx
GDICTH         = $(GDICT:.cxx=.h)
GDICTO         = $(patsubst %.cxx,tgt_$(ALICE_TARGET)/%.o,$(GDICT))

# Geant Objects

GOBJS          = $(patsubst %.cxx,tgt_$(ALICE_TARGET)/%.o,$(GSRCS)) $(FOBJS) $(GDICTO)

# FORTRAN Objectrs

FOBJS         = $(patsubst %.F,tgt_$(ALICE_TARGET)/%.o,$(FSRCS))

# C Objects

COBJS         = $(patsubst %.c,tgt_$(ALICE_TARGET)/%.o,$(CSRCS))

# C++ Objects

OBJS          = $(patsubst %.cxx,tgt_$(ALICE_TARGET)/%.o,$(SRCS)) $(DICTO)

# C++ Headers

DHDRS          = AliGeant3LinkDef.h

# Geant Objects

DOBJS          = $(patsubst %.cxx,tgt_$(ALICE_TARGET)/%.o,$(DSRCS)) $(DDICTO)

# C++ compilation flags

CXXFLAGS      = $(CXXOPTS) -I$(ROOTSYS)/include -I. -I$(ALICE_ROOT)/include/ \
	        -I$(G3INSTALL)/geant321 -I$(G3INSTALL)/geant3mc
 
# FORTRAN compilation flags

FFLAGS      = $(FOPT) -I$(ALICE_ROOT)/GEANT321

##### TARGETS #####
 
# Target

SLIBRARY	= $(LIBDIR)/libAliGeant3.$(SL) 
ALIBRARY	= $(LIBDIR)/libAliGeant3.a

default:	$(SLIBRARY)

$(LIBDIR)/libAliGeant3.$(SL):		$(GOBJS)

$(GDICT):				$(GHDRS)

$(DDICT):				$(DHDRS)

DICT	=		$(GDICT) $(DDICT)

DEPINC	+= -I$(ALICE_ROOT)/GEANT321

depend:			$(GSRCS) $(FSRCS)

TOCLEAN		= $(GOBJS) $(JOBJS) $(DOBJS) *Cint.cxx *Cint.h

CHECKS        = $(patsubst %.cxx,check/%.viol,$(GSRCS))

############################### General Macros ################################

include $(ALICE_ROOT)/conf/GeneralMacros

############################ Dependencies #####################################

-include tgt_$(ALICE_TARGET)/Make-depend 







