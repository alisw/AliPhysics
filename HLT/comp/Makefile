#
# Makefile for ROOT version of HLT code.
#
# Author: Anders Vestbo                       
#

CXX         = g++
LD          = g++
CXXFLAGS    = -O2 -Wall -fPIC -Woverloaded-virtual
LDFLAGS     = -O2
SOFLAGS     = -shared

TOPDIR = $(HOME)/level3code
OBJDIR = tgt_Linux

SRCS	= AliL3Modeller.cxx AliL3ModelTrack.cxx

DICT = AliL3CompCint.cxx
DICTH = AliL3CompCint.h
DICTO = $(OBJDIR)/AliL3CompCint.o

HDRS = $(SRCS:.cxx=.h) AliL3CompLinkDef.h

OBJS = $(patsubst %.cxx,$(OBJDIR)/%.o,$(SRCS)) $(DICTO)

LEVELSO = $(LEVEL3)/lib_$(USERNAME)/libAliL3Comp.so

all : $(OBJDIR) $(LEVELSO)


$(LEVELSO) : $(OBJS)
	$(LD) $(SOFLAGS) $(LDFLAGS) $^ -o $@

$(DICT) : $(HDRS)
	@echo "Generating dictionary..."
	@rootcint -f $(DICT) -c -I$(ROOTSYS) -I$(TOPDIR)/src/ $(HDRS)


$(OBJDIR)/%.o : %.cxx 
	$(CXX) $(CXXFLAGS) -c -I$(ROOTSYS)/include -I$(TOPDIR)/src -I$(TOPDIR)/hough -o $@ $<

$(OBJDIR) : 
	test ! -e $(OBJDIR) && mkdir -p $(OBJDIR)

clean :
	rm -f $(OBJDIR)/*.o
	rm -f $(LEVEL3)/lib_$(USERNAME)/libAliL3Comp.so
	rm $(DICT) $(DICTH) 
so :
	rm -f $(LEVEL3)/lib_$(USERNAME)/libAliL3Comp.so