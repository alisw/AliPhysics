
ROOTCINT	= rootcint
ROOTLIBS	= $(shell root-config --glibs) -lEG 
ROOTCFLAGS	= $(shell root-config --cflags)

#
# Compile variables 
#
LIBS		= $(ROOTLIBS) -L$(ALICE_ROOT)/lib/tgt_linux/ -lESD 
CPPFLAGS	= $(ROOTCFLAGS) -I$(ALICE_ROOT)/include -I$(ALICE_ROOT)/STEER/ -I./ 

CXX		= g++	
CXXFLAGS	= -c -g -Wall -fPIC 
LD		= g++
LDFLAGS		= -rdynamic -Wl,-rpath,./ $(LIBS) -o
SOFLAGS         = -shared -Wl,-soname,
CP		= cp

LIB	  = libESDtrackCuts.so
LIBO	  = AliESDtrackCutsCint.o \
	    AliESDtrackCuts.o 

#
# Rules

%Cint.cxx:%.h
	rootcint -f $@ -c $(ROOTCFLAGS) $(CPPFLAGS) -p -I./  $^

%.o:%.cxx
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $<

%.so:
	$(LD) $(SOFLAGS)$@ $^ -o $@

%:%.o
	$(LD) $(LDFLAGS) $(LIBS) $^ -o $@



#
#
#----------------------------------------------------------------------------


LIBRARIES = $(LIB)

#
# Targets 
#
.PHONY:all module 

all		: $(LIBRARIES) 

clean:
	@echo "Cleaning up ..."
	@$(RM) -f *~ core *Cint.* *.o

realclean: clean 
	@echo "Being very tidy ... "
	@$(RM) -f $(LIBRARIES)*
	@echo "done"


#
# Dependencies  
#
$(LIB):$(LIBO)



