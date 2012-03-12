PACKAGE = invalid-only-for-proof

include Makefile.arch

ifndef PACKCXXFLAGS
   PACKCXXFLAGS = $(CXXFLAGS)
endif

ALICEINC = -I.

### define include dir for local case and par case
ifneq ($(STEERBase_INCLUDE),)
  ALICEINC += -I../$(STEERBase_INCLUDE)
endif

ifneq ($(ESD_INCLUDE),)
  ALICEINC += -I../$(ESD_INCLUDE)
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

ifneq ($(EventMixing_INCLUDE),)
  ALICEINC += -I../$(EventMixing_INCLUDE)
endif

ifneq ($(JETAN_INCLUDE),)
  ALICEINC += -I../$(JETAN_INCLUDE)
endif

ifneq ($(JETANDEV_INCLUDE),)
  ALICEINC += -I../$(JETANDEV_INCLUDE)
endif

ifneq ($(CORRFW_INCLUDE),)
  ALICEINC += -I../$(CORRFW_INCLUDE)
endif

ifneq ($(PWGflowBase_INCLUDE),)
  ALICEINC += -I../$(PWGflowBase_INCLUDE)
endif

ifneq ($(PWGflowTasks_INCLUDE),)
  ALICEINC += -I../$(PWGflowTasks_INCLUDE)
endif

# only if no par file was loaded before
#ifeq ($(ALICEINC),-I.)
  ifneq ($(ALICE_ROOT),)
    ALICEINC += -I$(ALICE_ROOT)/include
    ALICEINC += -I$(ALICE_ROOT)/ANALYSIS
    ALICEINC += -I$(ALICE_ROOT)/PWG/FLOW/Base
    ALICEINC += -I$(ALICE_ROOT)/MUON
  endif
#endif

CXXFLAGS += $(ALICEINC) -g

SRCS = $(wildcard *.cxx)
HDRS = $(SRCS:.cxx=.h)
DHDR = $(PACKAGE)LinkDef.h

DHDR_PWG := $(DHDR)
HDRS_PWG := $(HDRS)
SRCS_PWG := $(SRCS) G__$(PACKAGE).cxx
OBJS_PWG := $(SRCS_PWG:.cxx=.o)

PARFILE       = $(PACKAGE).par

lib$(PACKAGE).so: $(OBJS_PWG)
	@echo "Linking" $@ ...
	@/bin/rm -f $@
	$(LD) $(SOFLAGS) $(LDFLAGS) $^ $(OutPutOpt) $@ $(EXPLLINKLIBS)
	@chmod a+x $@
	@echo "done"

%.o:    %.cxx %.h
	$(CXX) $(PACKCXXFLAGS) -c $< -o $@

clean:
	@rm -f $(OBJS_PWG) *.so G__$(PACKAGE).*

G__$(PACKAGE).cxx G__$(PACKAGE).h: $(HDRS) $(DHDR)
	@echo "Generating dictionaries ..." $(ALICEINC)
	rootcint -f $@ -c $(CINTFLAGS) $(ALICEINC) $^
