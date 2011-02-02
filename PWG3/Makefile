# $Id$

PACKAGE = PWG3

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
	ifneq ($(PWG3base_INCLUDE),)
		ALICEINC += -I../$(PWG3base_INCLUDE)
	endif
else
  ifneq ($(ALICE_ROOT),)
    ALICEINC += -I$(ALICE_ROOT)/include -I$(ALICE_ROOT)/PWG3/base -I$(ALICE_ROOT)/MUON -I$(ALICE_ROOT)/MUON/mapping
else
  endif
endif

SRCS     = $(wildcard $(subst PWG3,,$(PACKAGE))/*.cxx)
HDRS     = $(SRCS:.cxx=.h)
DHDR     = $(PACKAGE)LinkDef.h
CXXFLAGS     += $(ALICEINC) -g

DHDR_PWG3 := $(DHDR)
HDRS_PWG3 := $(HDRS)
SRCS_PWG3 := $(SRCS) G__$(PACKAGE).cxx
OBJS_PWG3 := $(SRCS_PWG3:.cxx=.o)

PARFILE       = $(PACKAGE).par

lib$(PACKAGE).so: $(OBJS_PWG3)
	@echo "Linking" $@ ...
	@/bin/rm -f $@
ifeq ($(ARCH),macosx)
	@$(LD) -bundle -undefined $(UNDEFOPT) $(LDFLAGS) $^ -o $@
else
	@$(LD) $(SOFLAGS) $(LDFLAGS) $^ -o $@
endif
	@chmod a+x $@
	@echo "done"

%.o:    %.cxx %.h
	$(CXX) $(CXXFLAGS) -c $< -o $@

G__$(PACKAGE).cxx G__$(PACKAGE).h: $(HDRS) $(DHDR)
	@echo "Generating dictionary ..."
	rootcint -f $@ -c $(ALICEINC) $^

clean:
	@rm -f $(OBJS)  *.so G__$(PACKAGE).*  $(PARFILE)
