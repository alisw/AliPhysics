############################### TPC Makefile ##################################

# Include machine specific definitions

include $(ALICE_ROOT)/conf/GeneralDef
include $(ALICE_ROOT)/conf/MachineDef.$(ALICE_TARGET)

PACKAGE = AliL3

# C++ sources


SRCS	       = AliL3ConfMapper.cxx AliL3ConfMapPoint.cxx AliLevel3.cxx\
	         AliL3ConfMapTrack.cxx AliL3ConfMapFit.cxx AliL3Track.cxx\
		 AliL3Transform.cxx AliL3TrackMerger.cxx AliL3VertexArray.cxx \
                 AliL3Vertex.cxx AliL3VertexFinder.cxx AliL3HoughTrack.cxx\
                 AliL3Merger.cxx AliL3GlobalMerger.cxx AliL3ClustFinder.cxx \
                 AliL3TrackArray.cxx AliL3InterMerger.cxx AliL3Logger.cxx \
                 AliL3MemHandler.cxx AliL3FileHandler.cxx AliL3Evaluate.cxx \
                 AliL3Benchmark.cxx AliL3Display.cxx
		 

# C++ Headers

HDRS          = $(SRCS:.cxx=.h) AliL3LinkDef.h
# Library dictionary

DICT          = AliL3Cint.cxx
DICTH         = $(DICT:.cxx=.h)
DICTO         = $(patsubst %.cxx,tgt_$(ALICE_TARGET)/%.o,$(DICT))

# FORTRAN Objectrs

FOBJS         = $(FSRCS:.f=.o)

# C Objects

COBJS         = $(CSRCS:.c=.o)

# C++ Objects

OBJS          = $(patsubst %.cxx,tgt_$(ALICE_TARGET)/%.o,$(SRCS)) $(DICTO)

# C++ compilation flags

CXXFLAGS      = $(CXXOPTS) -g -I$(ROOTSYS)/include -I. -I$(ALICE_ROOT)/include/  -I$(ALICE_ROOT)/TPC  -I$(ALICE_ROOT)/CONTAINERS

#CXXFLAGS      = $(CXXOPTS) -g -Wall -I$(ROOTSYS)/include -I. -I $(ALICE_ROOT)/TPC -I$(ALICE_ROOT)/include/ -DCOMPILING
 
# FORTRAN compilation flags

FFLAGS      = $(FOPT)
 
##### TARGETS #####
 
# Target

SLIBRARY	= $(ALICE)/mylibs/libAliL3.$(SL)

default:	$(SLIBRARY)

$(ALICE)/mylibs/libAliL3.$(SL):  $(OBJS)

$(DICT):		$(HDRS)

depend:			$(SRCS)

TOCLEAN			= $(OBJS) *Cint.h *Cint.cxx

############################### General Macros ################################

include $(ALICE_ROOT)/conf/GeneralMacros

############################ Dependencies #####################################

include tgt_$(ALICE_TARGET)/Make-depend 

###########
so:
	rm -fr $(ALICE)/mylibs/libAliL3.so
clean:
	rm -fr tgt_Linux/*.o
	rm -fr $(ALICE)/mylibs/libAliL3.so
	rm -fr $(DICT) $(DICTH) $(DICTO)
