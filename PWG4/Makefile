PACKAGE = invalid-only-for-proof

include $(ROOTSYS)/test/Makefile.arch
include lib$(PACKAGE).pkg

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

ifneq ($(JETAN_INCLUDE),)
  ALICEINC += -I../$(JETAN_INCLUDE)
endif

ifneq ($(CORRFW_INCLUDE),)
  ALICEINC += -I../$(CORRFW_INCLUDE)
endif

ifneq ($(PHOSUtils_INCLUDE),)
  ALICEINC += -I../$(PHOSUtils_INCLUDE)
  CXXFLAGS+=-D__PHOSUTIL__
endif

ifneq ($(EMCALUtils_INCLUDE),)
  ALICEINC += -I../$(EMCALUtils_INCLUDE)
  CXXFLAGS+=-D__EMCALUTIL__
endif

ifneq ($(PWG4PartCorrBase_INCLUDE),)
  ALICEINC += -I../$(PWG4PartCorrBase_INCLUDE)
endif

ifneq ($(PWG4PartCorrDep_INCLUDE),)
 ALICEINC += -I../$(PWG4PartCorrBase_INCLUDE) 
 ALICEINC += -I../$(PWG4PartCorrDep_INCLUDE)
endif

ifneq ($(PWG4GammaConv_INCLUDE),)
  ALICEINC += -I../$(PWG4GammaConv_INCLUDE)
endif

ifneq ($(PWG4JetTasks_INCLUDE),)
  ALICEINC += -I../$(PWG4JetTasks_INCLUDE)
endif

ifneq ($(JETAN_INCLUDE),)
  ALICEINC += -I../$(JETAN_INCLUDE)
endif

ifneq ($(PWG4JetCorrel_INCLUDE),)
  ALICEINC += -I../$(PWG4JetCorrel_INCLUDE)
endif

ifneq ($(PWG4JCORRAN_INCLUDE),)
  ALICEINC += -I../$(PWG4JCORRAN_INCLUDE)
endif

# only if no par file was loaded before
#ifeq ($(ALICEINC),-I.)
  ifneq ($(ALICE_ROOT),)
    ALICEINC += -I$(ALICE_ROOT)/include
    ALICEINC += -I$(ALICE_ROOT)/JETAN   # some extra includes
    ALICEINC += -I$(ALICE_ROOT)/PHOS
    ALICEINC += -I$(ALICE_ROOT)/EMCAL
  endif
#endif

CXXFLAGS += $(ALICEINC) -g

SRCS         += G__$(PACKAGE).cxx
OBJS          = $(SRCS:.cxx=.o)

lib$(PACKAGE).so: $(OBJS)
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
	$(CXX) $(PACKCXXFLAGS) -c $< -o $@

clean:
	@rm -f $(OBJS) *.so G__$(PACKAGE).*

G__$(PACKAGE).cxx G__$(PACKAGE).h: $(HDRS) $(DHDR)
	@echo "Generating dictionaries ..." $(ALICEINC)
	rootcint -f $@ -c $(CINTFLAGS) $(ALICEINC) $^

