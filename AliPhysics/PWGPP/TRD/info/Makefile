SrcSuf        = cxx
HdrSuf        = h
ExeSuf        =
ObjSuf        = o
DllSuf        = so
OutPutOpt     = -o

CXX           = g++
LD            = g++
CXXFLAGS      = -O -fPIC
LDFLAGS       = -O -fPIC
SOFLAGS       = -shared

ROOTCXX       = $(shell root-config --cflags)
ROOTLIBS      = $(shell root-config --libs)

ALICXX        = -I$(ALICE_ROOT)/include -I$(ALICE_ROOT)/TRD 

CXXFLAGS     += $(ROOTCXX) $(ALICXX)

# define module specific variables
FILE_LIST =  $(shell ls -1 ./*.$(SrcSuf))
FILES = $(basename $(FILE_LIST))
DICTIONARIES = trdTrackInfoDict
OBJECTS = $(addsuffix .$(ObjSuf),$(FILES))
OBJECTS += ./$(DICTIONARIES).$(ObjSuf)

#define headers
HDRS = $(addsuffix .$(HdrSuf),$(FILES)) 
HEADERS = $(notdir $(HDRS))

# define libs on which the main lib depends ! (this are defined in config/Makefile.flags)
LIBSDEPEND = $(ROOTLIBS) $(ALILIBS)
# define libs build by module
LIBS = libTRDtrackInfo.so


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

clean:
	@rm -fv $(DICTIONARIES)
	@rm -fv $(OBJECTS)
	@rm -fv $(LIBS)
	@rm -fv $(EXECS)
	@if [ "$(LIBS)" != "" ]; then rm -fv lib/$(LIBS); fi
	@if [ "$(EXECS)" != "" ]; then rm -fv bin/$(EXECS); fi
	

