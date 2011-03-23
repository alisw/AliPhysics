# AliRoot stat project
# Makefile for main library - libKDTree.so
#
# Authors
#       Marian Ivanov [M.Ivanov@gsi.de]
#       Alexandru Bercuci [A.Bercuci@gsi.de]
# Data     01.08.2007
# version  0.0 [pre release]
#

SrcSuf    = cxx
HdrSuf    = h
ExeSuf    =
ObjSuf    = o
DllSuf    = so
OutPutOpt = -o

CXX       = g++
LD        = g++
CXXFLAGS  = -O -g -pg -fPIC -Wall -Weffc++
LDFLAGS   = -O -g -pg -fPIC -Wall -Weffc++
SOFLAGS   = -shared -g -pg -fPIC

ROOTCXX  = $(shell root-config --cflags)
ROOTLIBS = $(shell root-config --libs) -lMinuit

MCCXX    = -I$(ALICE_ROOT)/STEER
MCLIBS   =

CXXFLAGS += $(ROOTCXX) $(MCCXX) $(LOCALCXX)

# define module specific variables
FILE_LIST =  $(shell ls -1 ./*.$(SrcSuf))
FILES = $(basename $(FILE_LIST))
DICTIONARIES = StatDict
OBJECTS = $(addsuffix .$(ObjSuf),$(FILES))
OBJECTS += ./$(DICTIONARIES).$(ObjSuf)

#define headers
HDRS = $(addsuffix .$(HdrSuf),$(FILES)) 
HEADERS = $(notdir $(HDRS))

# define libs on which the main lib depends ! (this are defined in config/Makefile.flags)
LIBSDEPEND = $(ROOTLIBS) $(MCLIBS)
# define libs build by module
LIBS = libSTAT.so

# rule for building executables
$(EXECS):	$(OBJECTS)
	@echo -e "\E[31mBuild executable: \E[1;31m$@\E[0m"
	@$(LD) $(LIBSDEPEND) $^ -o $@
	
# rule for building libraries
%.$(DllSuf):	$(OBJECTS)
	@echo -e "\E[31mBuild library: \E[1;31m$@\E[0m"
	@$(LD) $(SOFLAGS) $(LIBSDEPEND) $^ -o $@
	
# rule for building objects
%.$(ObjSuf):	%.$(SrcSuf)
	@echo -e "\E[31mCompile : \E[1;31m$^\E[0m"
	@$(CXX) $(CXXFLAGS) -c $< -o $@

#rule for building dictionary
%Dict.$(SrcSuf): %LinkDef.h
	@echo -e "\E[31mGenerate dictionary : \E[1;31m$@\E[0m"
	@rootcint -f $@ -c $(CXXFLAGS) $(HEADERS) $^
	
all: $(OBJECTS) $(LIBS) $(EXECS)

headers:
	@if [ "$(HEADERS)" != "" ]; then \
		if [ ! -d include ]; then mkdir include; fi; \
		cp -f $(HEADERS) include; \
	fi
	
clean:
	@rm -fv $(DICTIONARIES)
	@rm -fv $(OBJECTS)
	@find . -name "*~" -a -exec rm -fv \{} \;
	@if [ "$(LIBS)" != "" ]; then rm -fv ./$(LIBS); fi
	@if [ "$(EXECS)" != "" ]; then rm -fv ./$(EXECS); fi

ifeq ($(PLATFORM),macosx)
# We need to make both the .dylib and the .so
		$(LD) $(SOFLAGS)$@ $(LDFLAGS) $^ $(OutPutOpt) $@
ifneq ($(subst $(MACOSX_MINOR),,1234),1234)
ifeq ($(MACOSX_MINOR),4)
		ln -sf $@ $(subst .$(DllSuf),.so,$@)
else
		$(LD) -bundle -undefined $(UNDEFOPT) $(LDFLAGS) $^ \
		   $(OutPutOpt) $(subst .$(DllSuf),.so,$@)
endif
endif
else
ifeq ($(PLATFORM),win32)
		bindexplib $* $^ > $*.def
		lib -nologo -MACHINE:IX86 $^ -def:$*.def \
		   $(OutPutOpt)$(EVENTLIB)
		$(LD) $(SOFLAGS) $(LDFLAGS) $^ $*.exp $(LIBS) \
		   $(OutPutOpt)$@
else
		$(LD) $(SOFLAGS) $(LDFLAGS) $^ $(OutPutOpt) $@ $(EXPLLINKLIBS)
endif
endif
