PACKAGE = invalid-only-for-proof

include Makefile.arch

SRCS     = $(wildcard *.cxx)
HDRS     = $(SRCS:.cxx=.h)
DHDR     = $(PACKAGE)LinkDef.h


ifndef PACKCXXFLAGS
   PACKCXXFLAGS = $(CXXFLAGS)
endif

ALICEINC = -I.
ALICEINC += -I./$(PACKAGE)/

ifneq ($(ESD_INCLUDE),)
   ALICEINC += -I../$(ESD_INCLUDE)
endif

ifneq ($(AOD_INCLUDE),)
   ALICEINC += -I../$(AOD_INCLUDE)
endif

ifneq ($(STEERBase_INCLUDE),)
   ALICEINC += -I../$(STEERBase_INCLUDE)
endif

ifneq ($(ANALYSIS_INCLUDE),)
   ALICEINC += -I../$(ANALYSIS_INCLUDE)
endif

ifneq ($(OADB_INCLUDE),)
   ALICEINC += -I../$(OADB_INCLUDE)
endif

ifneq ($(ANALYSISalice_INCLUDE),)
   ALICEINC += -I../$(ANALYSISalice_INCLUDE)
endif

ifneq ($(CORRFW_INCLUDE),)
   ALICEINC += -I../$(CORRFW_INCLUDE)
endif

ifneq ($(PWGTools_INCLUDE),)
   ALICEINC += -I../$(PWGTools_INCLUDE)
endif

ifneq ($(PWGCFebye_INCLUDE),)
   ALICEINC += -I../$(PWGCFebye_INCLUDE) -I../$(PWGCFebye_INCLUDE)/LRC 
else
   ALICEINC += -I$(ALICE_ROOT)/PWGCF/EBYE -I$(ALICE_ROOT)/PWGCF/EBYE/LRC
endif

ifneq ($(PWGCFfemtoscopy_INCLUDE),)
   ALICEINC += -I../$(PWGCFfemtoscopy_INCLUDE) -I../$(PWGCFfemtoscopy_INCLUDE)/AliFemto
else
   ALICEINC += -I$(ALICE_ROOT)/PWGCF/FEMTOSCOPY -I$(ALICE_ROOT)/PWGCF/FEMTOSCOPY/AliFemto
endif

ifneq ($(PWGCFfemtoscopyUser_INCLUDE),)
   ALICEINC += -I../$(PWGCFfemtoscopyUser_INCLUDE)
endif

ifneq ($(PWGCFunicor_INCLUDE),)
   ALICEINC += -I../$(PWGCFunicor_INCLUDE)
endif

ifneq ($(PWGCFCorrelationsBase_INCLUDE),)
   ALICEINC += -I../$(PWGCFCorrelationsBase_INCLUDE)
endif

ifneq ($(PWGCFCorrelationsDPhi_INCLUDE),)
   ALICEINC += -I../$(PWGCFCorrelationsDPhi_INCLUDE)
endif

# only if no par file was loaded before
#ifeq ($(ALICEINC),-I.)
  ifneq ($(ALICE_ROOT),)
    ALICEINC += -I$(ALICE_ROOT)/include
  endif
#endif




CXXFLAGS += $(ALICEINC) -g

SRCS         += G__$(PACKAGE).cxx
OBJS          = $(SRCS:.cxx=.o)

lib$(PACKAGE).so: $(OBJS)
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

