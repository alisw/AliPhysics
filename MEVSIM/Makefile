############################### MEVSIM Makefile ###############################

# Include machine specific definitions

include $(ALICE_ROOT)/conf/GeneralDef
include $(ALICE_ROOT)/conf/MachineDef.$(ALICE_TARGET)

PACKAGE = MEVSIM


# C++ sources


##### MACROS #####

FSRCS	=  multiplicity_gen.F

FOBJS	= $(patsubst %.F,tgt_$(ALICE_TARGET)/%.o,$(FSRCS))

SRCS	= $(FSRCS) $(CSRCS)
OBJS	= $(FOBJS)

DSRCS	= dummymevsim.F
DOBJS	= $(patsubst %.F,tgt_$(ALICE_TARGET)/%.o,$(DSRCS))

# C++ compilation flags

CXXFLAGS      = $(CXXOPTS) $(CLIBCXXOPTS) $(CLIBDEFS)

 
# C compilation flags

CFLAGS      = $(COPT) $(CLIBCOPT) $(CLIBDEFS)
 
# FORTRAN compilation flags

FFLAGS      = $(FOPT) $(CLIBFOPT) $(CLIBDEFS)

##### TARGETS #####
 
# Target

SLIBRARY	= $(LIBDIR)/libmevsim.$(SL) $(LIBDIR)/libdummymevsim.$(SL)
ALIBRARY	= $(LIBDIR)/libmevsim.a

default:	$(SLIBRARY)

$(LIBDIR)/libmevsim.$(SL):  $(OBJS)
$(LIBDIR)/libdummymevsim.$(SL):  $(DOBJS)

depend:			$(CSRCS) $(DSRCS) $(SRCS)

TOCLEAN		= $(OBJS) $(DOBJS) *Cint.cxx *Cint.h

############################### General Macros ################################

include $(ALICE_ROOT)/conf/GeneralMacros

############################ Dependencies #####################################

-include tgt_$(ALICE_TARGET)/Make-depend 









