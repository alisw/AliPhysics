# $Id$

PACKAGE = PWGHF

include Makefile.arch

ALICEINC      = -I.

### define include dir for local case and par case
ifneq ($(ESD_INCLUDE),)
   ALICEINC += -I../$(ESD_INCLUDE)
	ifneq ($(STEERBase_INCLUDE),)
		ALICEINC += -I../$(STEERBase_INCLUDE)
	endif				
	ifneq ($(CDB_INCLUDE),)
		ALICEINC += -I../$(CDB_INCLUDE)
	endif				
	ifneq ($(AOD_INCLUDE),)
		ALICEINC += -I../$(AOD_INCLUDE)
	endif
	ifneq ($(ANALYSIS_INCLUDE),)
		ALICEINC += -I../$(ANALYSIS_INCLUDE)
	endif
	ifneq ($(ANALYSISalice_INCLUDE),)
		ALICEINC += -I../$(ANALYSISalice_INCLUDE)
	endif
	ifneq ($(CORRFW_INCLUDE),)
		ALICEINC += -I../$(CORRFW_INCLUDE)
	endif
	ifneq ($(PWGHFbase_INCLUDE),)
		ALICEINC += -I../$(PWGHFbase_INCLUDE)
	endif
else
  ifneq ($(ALICE_ROOT),)
    ALICEINC += -I$(ALICE_ROOT)/include -I$(ALICE_ROOT)/PWGHF/base
else
  endif
endif

SRCS     = $(wildcard $(subst PWGHF,,$(PACKAGE))/*.cxx)
HDRS     = $(SRCS:.cxx=.h)
DHDR     = $(PACKAGE)LinkDef.h
CXXFLAGS     += $(ALICEINC) -g

DHDR_PWGHF := $(DHDR)
HDRS_PWGHF := $(HDRS)
SRCS_PWGHF := $(SRCS) G__$(PACKAGE).cxx
OBJS_PWGHF := $(SRCS_PWGHF:.cxx=.o)

PARFILE       = $(PACKAGE).par

lib$(PACKAGE).so: $(OBJS_PWGHF)
	@echo "Linking" $@ ...
	@/bin/rm -f $@
	$(LD) $(SOFLAGS) $(LDFLAGS) $^ $(OutPutOpt) $@ $(EXPLLINKLIBS)
	@chmod a+x $@
	@echo "done"

%.o:    %.cxx %.h
	$(CXX) $(CXXFLAGS) -c $< -o $@

G__$(PACKAGE).cxx G__$(PACKAGE).h: $(HDRS) $(DHDR)
	@echo "Generating dictionary ..."
	rootcint -f $@ -c $(ALICEINC) $^

clean:
	@rm -f $(OBJS)  *.so G__$(PACKAGE).*  $(PARFILE)
