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

FOBJS	= $(FSRCS:.F=.o)

SRCS	= $(FSRCS)
OBJS	= $(FOBJS)

DSRCS	= dummypythia.F

DOBJS	= $(DSRCS:.F=.o)

# C++ compilation flags

CXXFLAGS      = $(CXXOPTS) $(CLIBCXXOPTS) $(CLIBDEFS)
 
# C compilation flags

CFLAGS      = $(COPT) (CLIBCOPT) $(CLIBDEFS)
 
# FORTRAN compilation flags

FFLAGS      = $(FOPT) $(CLIBFOPT) $(CLIBDEFS)
 
##### TARGETS #####
 
# Target

SLIBRARY	= $(ALICE_ROOT)/lib/libpythia.$(SL) $(ALICE_ROOT)/lib/libdummypythia.$(SL)
ALIBRARY	= $(ALICE_ROOT)/lib/libpythia.a

default:	$(SLIBRARY)

$(ALICE_ROOT)/lib/libpythia.$(SL):  $(OBJS)
$(ALICE_ROOT)/lib/libdummypythia.$(SL):  $(DOBJS)

depend:			$(SRCS) $(DSRCS)

TOCLEAN		= $(OBJS) $(DOBJS) *Cint.cxx *Cint.h

############################### General Macros ################################

include $(ALICE_ROOT)/conf/GeneralMacros

############################ Dependencies #####################################

include Make-depend 
