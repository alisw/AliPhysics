############################### TPC Makefile ##################################

# Include machine specific definitions

include $(ALICE_ROOT)/conf/GeneralDef
include $(ALICE_ROOT)/conf/MachineDef.$(ALICE_TARGET)

PACKAGE = AliL3Hough

# C++ sources


SRCS	       = AliL3HoughTransformer.cxx AliL3HoughPixel.cxx  \
		 AliL3HoughMaxFinder.cxx AliL3HoughEval.cxx AliL3HoughMerge.cxx 
		 

# C++ Headers

HDRS          = $(SRCS:.cxx=.h) AliL3Defs.h AliL3HoughLinkDef.h
# Library dictionary

DICT          = AliL3HoughCint.cxx
DICTH         = $(DICT:.cxx=.h)
DICTO         = $(patsubst %.cxx,tgt_$(ALICE_TARGET)/%.o,$(DICT))

# FORTRAN Objectrs

FOBJS         = $(FSRCS:.f=.o)

# C Objects

COBJS         = $(CSRCS:.c=.o)

# C++ Objects

OBJS          = $(patsubst %.cxx,tgt_$(ALICE_TARGET)/%.o,$(SRCS)) $(DICTO)

# C++ compilation flags

CXXFLAGS      = $(CXXOPTS) -g -I$(ROOTSYS)/include -I. -I$(ALICE_ROOT)/include/ -I$(ALICE_ROOT)/TPC \
		-I$(ALICE_ROOT)/CONTAINERS -I$(ALICE)/level3code/src 	

#CXXFLAGS      = $(CXXOPTS) -g -Wall -I$(ROOTSYS)/include -I. -I $(ALICE_ROOT)/TPC -I$(ALICE_ROOT)/include/ -DCOMPILING
 
# FORTRAN compilation flags

FFLAGS      = $(FOPT)
 
##### TARGETS #####
 
# Target

SLIBRARY	= $(ALICE)/mylibs/libAliL3Hough.$(SL)

default:	$(SLIBRARY)

$(ALICE)/mylibs/libAliL3Hough.$(SL):  $(OBJS)

$(DICT):		$(HDRS)

depend:			$(SRCS)

TOCLEAN			= $(OBJS) *Cint.h *Cint.cxx

############################### General Macros ################################

include $(ALICE_ROOT)/conf/GeneralMacros

############################ Dependencies #####################################

include tgt_$(ALICE_TARGET)/Make-depend 

###########
so:
	rm -fr $(ALICE)/mylibs/libAliL3Hough.so
clean:
	rm -fr tgt_Linux/*.o
	rm -fr $(ALICE)/mylibs/libAliL3Hough.so
	rm -fr $(DICT) $(DICTH) $(DICTO)
