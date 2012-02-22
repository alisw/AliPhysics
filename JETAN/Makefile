PACKAGE = invalid-only-for-proof

include Makefile.arch

SRCS     = $(wildcard *.cxx)
HDRS     = $(SRCS:.cxx=.h)
DHDR     = $(PACKAGE)LinkDef.h


ALICEINC      = -I.
ifneq ($(findstring dev,$(PACKAGE)),)
 ALICEINC     += -I./DEV
endif

### define include dir for local case and par case
ifneq ($(OADB_INCLUDE),)
   ALICEINC += -I../$(OADB_INCLUDE)
endif
ifneq ($(EMCALUtils_INCLUDE),)
   ALICEINC += -I../$(EMCALUtils_INCLUDE)
endif
ifneq ($(JETANDEV_INCLUDE),)
  ifneq ($(FASTJET_ROOT),)
   ALICEINC += -I$(FASTJET_ROOT)/include
  else
   ALICEINC += -I$(ALICE_ROOT)/JETAN/fastjet
  endif
endif
ifneq ($(JETAN_INCLUDE),)
 ifeq ($(findstring dev,$(PACKAGE)),)
  ifneq ($(FASTJET_ROOT),)
   ALICEINC += -I$(FASTJET_ROOT)/include
  else
   ALICEINC += -I$(ALICE_ROOT)/JETAN/fastjet
  endif
 endif
endif

ifneq ($(ESD_INCLUDE),)
   ALICEINC += -I../$(ESD_INCLUDE)
   ALICEINC += -I../$(AOD_INCLUDE)
   ALICEINC += -I../$(STEERBase_INCLUDE)
   ALICEINC += -I../$(ANALYSIS_INCLUDE)
   ALICEINC += -I../$(ANALYSISalice_INCLUDE)
   ifneq ($(JETANDEV_INCLUDE),)	
     ALICEINC += -I../$(JETANDEV_INCLUDE)
   endif
   ifneq ($(JETAN_INCLUDE),)
    ifeq ($(findstring dev,$(PACKAGE)),)
     ALICEINC += -I../$(JETAN_INCLUDE)
    endif
   endif
 else
    ifneq ($(ALICE_ROOT),)
      ifneq ($(findstring dev,$(PACKAGE)),)
       ALICEINC += -I$(ALICE_ROOT)/include -I$(ALICE_ROOT)/JETAN/DEV -I$(ALICE_ROOT)/JETAN/fastjet -I$(ALICE_ROOT)/ANALYSIS
      else
       ALICEINC += -I$(ALICE_ROOT)/include -I$(ALICE_ROOT)/JETAN -I$(ALICE_ROOT)/JETAN/fastjet -I$(ALICE_ROOT)/ANALYSIS
     endif
    endif
endif

CXXFLAGS     += $(ALICEINC) -g


SRCS         += G__$(PACKAGE).cxx
OBJS          = $(SRCS:.cxx=.o)

lib$(PACKAGE).so: $(OBJS)
	@echo "Linking" $@ ...
	@/bin/rm -f $@
	$(LD) $(SOFLAGS) $(LDFLAGS) $^ $(OutPutOpt) $@ $(EXPLLINKLIBS)
	@chmod a+x $@
	@echo "done"

%.o:    %.cxx %.h
	$(CXX) $(CXXFLAGS) -c $< -o $@

clean:
	@rm -f $(OBJS) *.so G__$(PACKAGE).*

G__$(PACKAGE).cxx G__$(PACKAGE).h: $(HDRS) $(DHDR)
	@echo "Generating dictionary ..."
	rootcint -f $@ -c $(CINTFLAGS) $(ALICEINC) $^

