#
# Config variables
#

MAJOR		= 0
MINOR		= 1

#PROGRAM		= maketracktree
#PROGOBJS	= maketracktree.o
PROGRAM		= 
PROGOBJS	= 

PROGRAM2	= 
PROGOBJS2	= 

LIBNAME		= MyDeDxAnalysis
LIBRARY		= lib$(LIBNAME).so
LIBOBJS		= AliHighPtDeDxBase.o	\
		  AliHighPtDeDxCalib.o	\
		  AliHighPtDeDxData.o	\
		  AliHighPtDeDxMc.o	\
		  AliHighPtDeDxSpectra.o	\
		  AliHighPtDeDxBaseCint.o


#
# ROOT and BRAT variables
#
ROOTCINT	= rootcint 
ROOTCFLAGS	= $(shell root-config --cflags)
#ROOTLIBS	= $(shell root-config --glibs) -lEG -lGeom -lTreePlayer
ROOTLIBS	= $(shell root-config --glibs)
ROOTINCS	= -I./ 
ACCPATH         = $(shell pwd)

#
#
#
ALILIBS         = -L$(ALICE_ROOT)/lib/tgt_linux/ -lESD
#LIBS		= -L./ -l$(LIBNAME) $(ROOTLIBS) $(ALILIBS) -L$(HOME)/lib -lTOption
LIBS		= -L./ -l$(LIBNAME) $(ROOTLIBS)
ALIINCLUDE      = -I$(ALICE_ROOT)/include/ -I$(ALICE_ROOT)/STEER/ -I$(ALICE_ROOT)/TPC/
#INCLUDES	= -I./ $(ALIINCLUDE) -I$(HOME)/include/ 
INCLUDES	= -I./ 

# 
# Preprocessor, Compiler, and Linker variables 
#
CPPFLAGS	= $(INCLUDES) $(ROOTCFLAGS)
CXX		= g++
#CXXFLAGS	= -c -g -Wall -o  
CXXFLAGS	= -c -g -Wall -fPIC -o  
LD		= g++
LDFLAGS		= -rdynamic -Wl,-rpath,$(ACCPATH) $(LIBS) -o 
SOFLAGS         = -shared  -Wl,-soname,
LN		= ln -fs 
OUTOPT		= -o

#
# Rules
#
#%Cint.cxx:%Inc.h
%Cint.cxx:Include.h LinkDef.h
	rootcint -f $@ -c $(ROOTCFLAGS) -p -I./ $^ 
#$*LinkDef.h

%.o:%.cc
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $@ $<

%.o:%.cxx
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $@ $<

%.so:
	$(LD) $(SOFLAGS)$@ $(OUTOPT) $@ $^

%:%.o
	$(LD) $(LDFLAGS) $@ $^


#
# Targets 
#
all:	$(LIBRARY) $(PROGRAM) $(PROGRAM2) 

install: all
	@echo "cannot install, don't know your install dir"
clean:
	rm -f core *~ *.o *Cint.* $(LIBRARY)* $(PROGRAM) $(PROGRAM2) 

realclean: clean
	rm -f $(LIBRARY)* $(PROGRAM) $(PROGRAM2) 

#
# Dependencies 
#
$(LIBRARY):$(LIBOBJS)
$(PROGRAM):$(PROGOBJS)
$(PROGRAM2):$(PROGOBJS2)


