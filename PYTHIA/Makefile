############################### PYTHIA Makefile ###############################

# Include machine specific definitions

include $(ALICE_ROOT)/conf/GeneralDef
include $(ALICE_ROOT)/conf/MachineDef.$(ALICE_TARGET)

PACKAGE = PYTHIA

# C++ sources

SRCS          =
##### MACROS #####

FSRCS	= $(filter-out %pdfset.F,$(filter-out %structm.F,$(wildcard pythia/*.F))) \
          $(wildcard jetset/*.F)

FOBJS	= $(patsubst %.F,tgt_$(ALICE_TARGET)/%.o,$(FSRCS))

SRCS	= $(FSRCS)
OBJS	= $(FOBJS)

DSRCS	= dummypythia.F

DOBJS	= $(patsubst %.F,tgt_$(ALICE_TARGET)/%.o,$(DSRCS))

# C++ compilation flags

CXXFLAGS      = $(CXXOPTS) $(CLIBCXXOPTS) $(CLIBDEFS)
 
# C compilation flags

CFLAGS      = $(COPT) (CLIBCOPT) $(CLIBDEFS)
 
# FORTRAN compilation flags

FFLAGS      = $(FOPT) $(CLIBFOPT) $(CLIBDEFS)
 
##### TARGETS #####
 
# Target

SLIBRARY	= $(LIBDIR)/libpythia.$(SL) $(LIBDIR)/libdummypythia.$(SL)
ALIBRARY	= $(LIBDIR)/libpythia.a

default:	$(SLIBRARY)

$(LIBDIR)/libpythia.$(SL):  $(OBJS)
$(LIBDIR)/libdummypythia.$(SL):  $(DOBJS)

depend:			$(SRCS) $(DSRCS)

TOCLEAN		= $(OBJS) $(DOBJS) *Cint.cxx *Cint.h

############################### General Macros ################################

include $(ALICE_ROOT)/conf/GeneralMacros

############################ Dependencies #####################################

-include Make-depend 
