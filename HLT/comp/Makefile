#$Id$
###########################
# Makefile for HLT code.  #
#                         #
# Author: Anders Vestbo   #                    
###########################

CXX         = g++
LD          = g++
CXXFLAGS    = -O2 -Wall -fPIC -Woverloaded-virtual
LDFLAGS     = -O2
SOFLAGS     = -shared

TOPDIR = $(HOME)/level3code
OBJDIR = lib
LIBDIR = $(LEVEL3)/lib_$(USERNAME)

SRCS	= AliL3Modeller.cxx AliL3ModelTrack.cxx AliL3Compress.cxx

DICT = AliL3CompCint.cxx
DICTH = AliL3CompCint.h
DICTO = $(OBJDIR)/AliL3CompCint.o

HDRS = $(SRCS:.cxx=.h) bitio.h errhand.h AliL3CompLinkDef.h

OBJS = $(patsubst %.cxx,$(OBJDIR)/%.o,$(SRCS)) $(OBJDIR)/bitio.o $(OBJDIR)/errhand.o $(DICTO)

LEVELSO = $(LIBDIR)/libAliL3Comp.so

all : $(OBJDIR) $(LIBDIR) $(LEVELSO)


$(LEVELSO) : $(OBJS)
	$(LD) $(SOFLAGS) $(LDFLAGS) $^ -o $@

$(DICT) : $(HDRS)
	@echo "Generating dictionary..."
	@rootcint -f $(DICT) -c -I$(ROOTSYS) -I$(TOPDIR)/src/ $(HDRS)


$(OBJDIR)/%.o : %.cxx 
	$(CXX) $(CXXFLAGS) -c -I$(ROOTSYS)/include -I$(TOPDIR)/src -I$(TOPDIR)/hough -o $@ $<

$(OBJDIR) : 
	test ! -e $(OBJDIR) && mkdir -p $(OBJDIR)

$(LIBDIR) :
	test ! -e $(LIBDIR) && mkdir -p $(LIBDIR)

$(OBJDIR)/bitio.o : bitio.c
	$(CXX) -c -o $(OBJDIR)/bitio.o bitio.c
$(OBJDIR)/errhand.o : errhand.c
	$(CXX) -c -o $(OBJDIR)/errhand.o errhand.c

clean :
	rm -f $(OBJDIR)/*.o
	rm -f $(LIBDIR)/libAliL3Comp.so
	rm $(DICT) $(DICTH) 
so :
	rm -f $(LIBDIR)/libAliL3Comp.so