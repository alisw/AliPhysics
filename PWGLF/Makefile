PACKAGE = PWGLFresonances

include Makefile.arch

SRCS = 
HDRS = 
DHDR = 

ifndef PACKCXXFLAGS
   PACKCXXFLAGS = $(CXXFLAGS)
endif

PARPATH 	:= $(shell pwd)
ALICEINC 	= -I.
ALICEINC 	+= -I$(PACKAGE) -I$(PACKAGE)/RESONANCES
ALIPARLIBS 	= 

ifneq ($(STEERBase_INCLUDE),)
   ALICEINC += -I../$(STEERBase_INCLUDE)
   ALIPARLIBS += $(STEERBase_LIBS)
else
   ifneq ($(PACKAGE),STEERBase)
     ALICELIBS += -lSTEERBase
   else
     ALIPARLIBS += -L$(PARPATH) -lSTEERBase
   endif
endif

ifneq ($(ESD_INCLUDE),)
   ALICEINC += -I../$(ESD_INCLUDE)
   ALIPARLIBS += $(ESD_LIBS)
else
   ifneq ($(PACKAGE),ESD)
     ALICELIBS += -lESD
   else
     ALIPARLIBS += -L$(PARPATH) -lESD
   endif
endif

ifneq ($(AOD_INCLUDE),)
   ALICEINC += -I../$(AOD_INCLUDE)
   ALIPARLIBS += $(AOD_LIBS)
else
   ifneq ($(PACKAGE),AOD)
     ALICELIBS += -lAOD
   else
     ALIPARLIBS += -L$(PARPATH) -lAOD
   endif
endif

ifneq ($(ANALYSIS_INCLUDE),)
   ALICEINC += -I../$(ANALYSIS_INCLUDE)
   ALIPARLIBS += $(ANALYSIS_LIBS)
else
   ifneq ($(PACKAGE),ANALYSIS)
     ALICELIBS += -lANALYSIS
   else
     ALIPARLIBS += -L$(PARPATH) -lANALYSIS
   endif
endif

ifneq ($(OADB_INCLUDE),)
   ALICEINC += -I../$(OADB_INCLUDE)
   ALIPARLIBS += $(OADB_LIBS)
else
   ifneq ($(PACKAGE),OADB)
     ALICELIBS += -lOADB
   else
     ALIPARLIBS += -L$(PARPATH) -lOADB
   endif
endif

ifneq ($(ANALYSISalice_INCLUDE),)
   ALICEINC += -I../$(ANALYSISalice_INCLUDE)
   ALIPARLIBS += $(ANALYSISalice_LIBS)
else
   ifneq ($(PACKAGE),ANALYSISalice)
     ALICELIBS += -lANALYSISalice
   else
     ALIPARLIBS += -L$(PARPATH) -lANALYSISalice
   endif
endif

ifneq ($(EventMixing_INCLUDE),)
   ALICEINC += -I../$(EventMixng_INCLUDE)
   ALIPARLIBS += $(EventMixing_LIBS)
else
   ifneq ($(PACKAGE),EventMixing)
     ALICELIBS += -lEventMixing
   else
     ALIPARLIBS += -L$(PARPATH) -lEventMixing
   endif
endif

ifneq ($(CORRFW_INCLUDE),)
   ALICEINC += -I../$(CORRFW_INCLUDE)
   ALIPARLIBS += $(CORRFW_LIBS)
else
   ifneq ($(PACKAGE),CORRFW)
     ALICELIBS += -lCORRFW
   else
     ALIPARLIBS += -L$(PARPATH) -lCORRFW
   endif
endif

ifneq ($(PWGLFresonances_INCLUDE),)
   ALICEINC += -I../$(PWGLFresonances_INCLUDE)
   ALIPARLIBS += $(PWGLFresonances_LIBS)
else
   ifneq ($(PACKAGE),PWGLFresonances)
     ALICELIBS += -lPWGLFresonances
   else     
     ALIPARLIBS += -L$(PARPATH) -lPWGLFresonances
   endif
endif

# only if no par file was loaded before
ifeq ($(ALICEINC),-I.)
  ifneq ($(ALICE_ROOT),)
    ALICEINC += -I$(ALICE_ROOT)/include
    ALIPARLIBS += -L$(ALICE_ROOT)/lib $(ALICELIBS)
  endif
  ifneq ($(ALICE_PHYSICS),)
    ALICEINC += -I$(ALICE_PHYSICS)/include
    ALIPARLIBS += -L$(ALICE_PHYSICS)/lib $(ALICELIBS)
  endif
endif

ALICEINC 	+= -I$(ALICE_ROOT)/include -I${ALICE_PHYSICS}/include
ALIPARLIBS 	+= -L$(ALICE_ROOT)/lib $(ALICELIBS) \
	      	   -L${ALICE_PHYSICS}/lib 

ROOTLIBDIR    	:= $(shell $(RC) --libdir)
PACKLDFLAGS 	:= $(ALIPARLIBS) -L$(ROOTLIBDIR) \
		   -lMinuit -lTree -lXMLParser -lVMC -lMinuit -lEG -lGeom
CXXFLAGS 	+= $(ALICEINC) -g

SRCS         	+= G__$(PACKAGE).cxx
OBJS          	=  $(SRCS:.cxx=.o)

lib$(PACKAGE).so: $(OBJS)
	@echo "Linking" $@ ...
	@/bin/rm -f $@
ifeq ($(PLATFORM),macosx)
	$(LD) $(SOFLAGS) -flat_namespace -undefined suppress $(LDFLAGS) \
		$^ $(OutPutOpt) $@ $(PACKLDFLAGS) $(EXPLLINKLIBS)
else
	$(LD) $(SOFLAGS) $(LDFLAGS) $^ $(OutPutOpt) $@ \
		$(PACKLDFLAGS) $(EXPLLINKLIBS)
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

