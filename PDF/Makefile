############################### PDF Makefile ##################################

# Include machine specific definitions

include $(ALICE_ROOT)/conf/GeneralDef
include $(ALICE_ROOT)/conf/MachineDef.$(ALICE_TARGET)

PACKAGE = PDF

# C++ sources

SRCS          =

FSRCS	= $(wildcard spdf/*.F)

FOBJS	= $(patsubst %.F,tgt_$(ALICE_TARGET)/%.o,$(FSRCS))

SRCS	= $(FSRCS)
OBJS	= $(FOBJS)

# C++ compilation flags

CXXFLAGS      = $(CXXOPTS) $(CLIBCXXOPTS) $(CLIBDEFS) -I.
 
# C compilation flags

CFLAGS      = $(COPT) $(CLIBCOPT) $(CLIBDEFS) -I.
 
# FORTRAN compilation flags

FFLAGS      = $(FOPT) $(CLIBFOPT) $(CLIBDEFS) -I.
 
# Target

SLIBRARY	= $(LIBDIR)/libpdf.$(SL)
ALIBRARY	= $(LIBDIR)/libpdf.a

default:	$(SLIBRARY)

$(LIBDIR)/libpdf.$(SL):  $(OBJS)

depend:			$(SRCS)

TOCLEAN		= $(OBJS)

############################### General Macros ################################

include $(ALICE_ROOT)/conf/GeneralMacros

############################ Dependencies #####################################

-include Make-depend 
