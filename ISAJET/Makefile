############################### PDF Makefile ##################################

# Include machine specific definitions

include $(ALICE_ROOT)/conf/GeneralDef
include $(ALICE_ROOT)/conf/MachineDef.$(ALICE_TARGET)

PACKAGE = ISAJET

PDFDIRS = code isadata isarun isasusy isatape test openfile pdfinit  utils/cern_lib

# C++ sources

SRCS    =

FSRCS	= $(wildcard $(patsubst %,%/*.F,$(PDFDIRS)))

FOBJS	=  $(patsubst %.F,tgt_$(ALICE_TARGET)/%.o,$(FSRCS))

SRCS	= $(FSRCS)
OBJS	= $(FOBJS)

# C++ compilation flags

CXXFLAGS      = $(CXXOPTS) $(CLIBCXXOPTS) $(CLIBDEFS) -I.
 
# C compilation flags

CFLAGS      = $(COPT) $(CLIBCOPT) $(CLIBDEFS) -I.
 
# FORTRAN compilation flags

FFLAGS      = $(FOPT) $(CLIBFOPT) $(CLIBDEFS) -DCERNLIB_DOUBLE  -DCERNLIB_PDFLIB -I.
 
# Target

SLIBRARY	= $(LIBDIR)/libisajet.$(SL)
ALIBRARY	= $(LIBDIR)/libisajet.a

default:	$(SLIBRARY)

$(LIBDIR)/libisajet.$(SL):  $(OBJS)

depend:		$(SRCS)

TOCLEAN		= $(OBJS)

############################### General Macros ################################

include $(ALICE_ROOT)/conf/GeneralMacros

############################ Dependencies #####################################

-include tgt_$(ALICE_TARGET)/Make-depend 











