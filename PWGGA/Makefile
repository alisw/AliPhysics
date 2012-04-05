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

ifneq ($(PHOSUtils_INCLUDE),)
  ALICEINC += -I../$(PHOSUtils_INCLUDE)
  CXXFLAGS+=-D__PHOSUTIL__
endif

ifneq ($(EMCALUtils_INCLUDE),)
  ALICEINC += -I../$(EMCALUtils_INCLUDE)
  CXXFLAGS+=-D__EMCALUTIL__
endif

ifneq ($(PWGGACaloTasks_INCLUDE),)
  ALICEINC += -I../$(PWGGACaloTasks_INCLUDE)
endif

ifneq ($(PWGCaloTrackCorrBase_INCLUDE),)
  ALICEINC += -I../$(PWGCaloTrackCorrBase_INCLUDE)
endif

ifneq ($(PWGGACaloTrackCorrelations_INCLUDE),)
  ALICEINC += -I../$(PWGCaloTrackCorrBase_INCLUDE)
  ALICEINC += -I../$(PWGGACaloTrackCorrelations_INCLUDE)
endif

ifneq ($(PWGGAEMCALTasks_INCLUDE),)
  ALICEINC += -I../$(PWGGAEMCALTasks_INCLUDE)
endif

ifneq ($(PWGGAGammaConv_INCLUDE),)
  ALICEINC += -I../$(PWGGAGammaConv_INCLUDE)
endif

ifneq ($(PWGGAPHOSTasks_INCLUDE),)
  ALICEINC += -I../$(PWGGAPHOSTasks_INCLUDE)
endif

# only if no par file was loaded before
#ifeq ($(ALICEINC),-I.)
  ifneq ($(ALICE_ROOT),)
    ALICEINC += -I$(ALICE_ROOT)/include
    ALICEINC += -I$(ALICE_ROOT)/ANALYSIS
    ALICEINC += -I$(ALICE_ROOT)/PWG
    ALICEINC += -I$(ALICE_ROOT)/PHOS
    ALICEINC += -I$(ALICE_ROOT)/EMCAL
  endif
#endif

CXXFLAGS += $(ALICEINC) -g

SRCS     = $(wildcard $(subst PWGGA,,$(PACKAGE))/*.cxx)
HDRS     = $(SRCS:.cxx=.h)
DHDR     = $(PACKAGE)LinkDef.h

DHDR_PWGGA := $(DHDR)
HDRS_PWGGA := $(HDRS)
SRCS_PWGGA := $(SRCS) G__$(PACKAGE).cxx
OBJS_PWGGA := $(SRCS_PWGGA:.cxx=.o)

PARFILE       = $(PACKAGE).par

lib$(PACKAGE).so: $(OBJS_PWGGA)
	@echo "Linking" $@ ...
	@/bin/rm -f $@
	$(LD) $(SOFLAGS) $(LDFLAGS) $^ $(OutPutOpt) $@ $(EXPLLINKLIBS)
	@chmod a+x $@
	@echo "done"

%.o:    %.cxx %.h
	$(CXX) $(PACKCXXFLAGS) -c $< -o $@

clean:
	@rm -f $(OBJS) *.so G__$(PACKAGE).*

G__$(PACKAGE).cxx G__$(PACKAGE).h: $(HDRS) $(DHDR)
	@echo "Generating dictionaries ..." $(ALICEINC)
	rootcint -f $@ -c $(CINTFLAGS) $(ALICEINC) $^
