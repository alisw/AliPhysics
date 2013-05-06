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

ifneq ($(CORRFW_INCLUDE),)
  ALICEINC += -I../$(CORRFW_INCLUDE)
endif

ifneq ($(PWGTools_INCLUDE),)
  ALICEINC += -I../$(PWGTools_INCLUDE)
endif

ifneq ($(PHOSUtils_INCLUDE),)
  ALICEINC += -I../$(PHOSUtils_INCLUDE)
  CXXFLAGS+=-D__PHOSUTIL__
endif

ifneq ($(EMCALUtils_INCLUDE),)
  ALICEINC += -I../$(EMCALUtils_INCLUDE)
  CXXFLAGS+=-D__EMCALUTIL__
endif

ifneq ($(EMCALbase_INCLUDE),)
  ALICEINC += -I../$(EMCALbase_INCLUDE)
endif

ifneq ($(EMCALsim_INCLUDE),)
  ALICEINC += -I../$(EMCALsim_INCLUDE)
endif

ifneq ($(EMCALrec_INCLUDE),)
  ALICEINC += -I../$(EMCALrec_INCLUDE)
endif

ifneq ($(PWGCaloTrackCorrBase_INCLUDE),)
  ALICEINC += -I../$(PWGCaloTrackCorrBase_INCLUDE)
endif

ifneq ($(PWG4PartCorrDep_INCLUDE),)
 ALICEINC += -I../$(PWGCaloTrackCorrBase_INCLUDE)
 ALICEINC += -I../$(PWGGACaloTrackCorrelations_INCLUDE)
endif

ifneq ($(PWGGACaloTrackCorrelations_INCLUDE),)
  ALICEINC += -I../$(PWG4GammaConv_INCLUDE)
endif

ifneq ($(JETAN_INCLUDE),)
  ALICEINC += -I../$(JETAN_INCLUDE)
endif

ifneq ($(JETANDEV_INCLUDE),)
  ALICEINC += -I../$(JETANDEV_INCLUDE)
endif

ifneq ($(PWGEMCAL_INCLUDE),)
  ALICEINC += -I../$(PWGEMCAL_INCLUDE)
endif

ifneq ($(PWGJEEMCALJetTasks_INCLUDE),)
  ALICEINC += -I../$(PWGJEEMCALJetTasks_INCLUDE)
endif

# only if no par file was loaded before
#ifeq ($(ALICEINC),-I.)
  ifneq ($(ALICE_ROOT),)
    ALICEINC += -I$(ALICE_ROOT)/include
    ALICEINC += -I$(ALICE_ROOT)/ANALYSIS
    ALICEINC += -I$(ALICE_ROOT)/PHOS
    ALICEINC += -I$(ALICE_ROOT)/EMCAL
    ALICEINC += -I$(ALICE_ROOT)/JETAN
    ALICEINC += -I$(ALICE_ROOT)/JETAN/fastjet
    ALICEINC += -I$(ALICE_ROOT)/PWGJE/EMCALJetTasks
    ALICEINC += -I$(ALICE_ROOT)/PWG/CaloTrackCorrBase
    ALICEINC += -I$(ALICE_ROOT)/PWGHF/vertexingHF
  endif
#endif

CXXFLAGS += $(ALICEINC) -g

SRCS = $(wildcard *.cxx)
HDRS = $(SRCS:.cxx=.h)
DHDR = $(PACKAGE)LinkDef.h

DHDR_PWGJE := $(DHDR)
HDRS_PWGJE := $(HDRS)
SRCS_PWGJE := $(SRCS) G__$(PACKAGE).cxx
OBJS_PWGJE := $(SRCS_PWGJE:.cxx=.o)

PARFILE       = $(PACKAGE).par

lib$(PACKAGE).so: $(OBJS_PWGJE)
	@echo "Linking" $@ ...
	@/bin/rm -f $@
	$(LD) $(SOFLAGS) $(LDFLAGS) $^ $(OutPutOpt) $@ $(EXPLLINKLIBS)
	@chmod a+x $@
	@echo "done"

%.o:    %.cxx %.h
	$(CXX) $(PACKCXXFLAGS) -c $< -o $@

clean:
	@rm -f $(OBJS_PWGJE) *.so G__$(PACKAGE).*

G__$(PACKAGE).cxx G__$(PACKAGE).h: $(HDRS) $(DHDR)
	@echo "Generating dictionaries ..." $(ALICEINC)
	rootcint -f $@ -c $(CINTFLAGS) $(ALICEINC) $^
