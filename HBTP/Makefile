############################### HBTPROCESSOR Makefile ###############################

# Include machine specific definitions

include $(ALICE_ROOT)/conf/GeneralDef
include $(ALICE_ROOT)/conf/MachineDef.$(ALICE_TARGET)

PACKAGE = HBTPROCESSOR
MAKEFLAGS =

# C sources

CSRCS         = 


# C Objects

COBJS         = $(patsubst %.c,tgt_$(ALICE_TARGET)/%.o,$(CSRCS))

##### MACROS #####

FSRCS	=  hbt_event_processor.f
#hbt_event_processor.f

FOBJS	= $(patsubst %.f,tgt_$(ALICE_TARGET)/%.o,$(FSRCS))

SRCS	= $(FSRCS) $(CSRCS)
OBJS	= $(FOBJS) $(COBJS)

DSRCS	= dummyhbtproc.f
#dummyhbtprocessor.f

DOBJS	= $(patsubst %.f,tgt_$(ALICE_TARGET)/%.o,$(DSRCS))

# C++ compilation flags

CXXFLAGS      = $(CXXOPTS) $(CLIBCXXOPTS) $(CLIBDEFS)
 
# C compilation flags

CFLAGS      = $(COPT) $(CLIBCOPT) $(CLIBDEFS)
 
# FORTRAN compilation flags

FFLAGS      = $(FOPT) $(CLIBFOPT) $(CLIBDEFS)
 
##### TARGETS #####
 
# Target

SLIBRARY	= $(LIBDIR)/libHBTP.$(SL) $(LIBDIR)/libdummyHBTP.$(SL)
ALIBRARY	= $(LIBDIR)/libHBTP.a $(LIBDIR)/libdummyHBTP.a

default:	$(SLIBRARY)

$(LIBDIR)/libHBTP.$(SL):  $(OBJS)
$(LIBDIR)/libdummyHBTP.$(SL):  $(DOBJS)

depend:			$(CSRCS) $(DSRCS)

TOCLEAN		= $(OBJS) $(DOBJS) *Cint.cxx *Cint.h

############################### General Macros ################################

include $(ALICE_ROOT)/conf/GeneralMacros

############################ Dependencies #####################################

-include tgt_$(ALICE_TARGET)/Make-depend 

clean:
	rm -f $(SLIBRARY) $(OBJS) $(DOBJS)

