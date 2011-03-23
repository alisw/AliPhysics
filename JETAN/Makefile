PACKAGE = invalid-only-for-proof

include $(ROOTSYS)/test/Makefile.arch

SRCS     = $(wildcard *.cxx)
HDRS     = $(SRCS:.cxx=.h)
DHDR     = $(PACKAGE)LinkDef.h


ALICEINC      = -I.

### define include dir for local case and par case
ifneq ($(OADB_INCLUDE),)
   ALICEINC += -I../$(OADB_INCLUDE)
endif
ifneq ($(EMCALUtils_INCLUDE),)
   ALICEINC += -I../$(EMCALUtils_INCLUDE)
endif
ifneq ($(JETAN_INCLUDE),)
  ifneq ($(FASTJET_ROOT),)
   ALICEINC += -I$(FASTJET_ROOT)/include
  else
   ALICEINC += -I$(ALICE_ROOT)/JETAN/fastjet
  endif
endif
ifneq ($(ESD_INCLUDE),)
   ALICEINC += -I../$(ESD_INCLUDE)
   ALICEINC += -I../$(AOD_INCLUDE)
   ALICEINC += -I../$(STEERBase_INCLUDE)
   ALICEINC += -I../$(ANALYSIS_INCLUDE)
   ALICEINC += -I../$(ANALYSISalice_INCLUDE)
   ifneq ($(JETAN_INCLUDE),)	
     ALICEINC += -I../$(JETAN_INCLUDE)
   endif
 else
    ifneq ($(ALICE_ROOT),)
       ALICEINC += -I$(ALICE_ROOT)/include -I$(ALICE_ROOT)/JETAN -I$(ALICE_ROOT)/ANALYSIS
     endif
endif

CXXFLAGS     += $(ALICEINC) -g


SRCS         += G__$(PACKAGE).cxx
OBJS          = $(SRCS:.cxx=.o)

lib$(PACKAGE).so: $(OBJS)
	@echo "Linking" $@ ...
	@/bin/rm -f $@
ifeq ($(PLATFORM),macosx)
# We need to make both the .dylib and the .so
		$(LD) $(SOFLAGS)$@ $(LDFLAGS) $^ $(OutPutOpt) $@
ifneq ($(subst $(MACOSX_MINOR),,1234),1234)
ifeq ($(MACOSX_MINOR),4)
		ln -sf $@ $(subst .$(DllSuf),.so,$@)
else
		$(LD) -bundle -undefined $(UNDEFOPT) $(LDFLAGS) $^ \
		   $(OutPutOpt) $(subst .$(DllSuf),.so,$@)
endif
endif
else
ifeq ($(PLATFORM),win32)
		bindexplib $* $^ > $*.def
		lib -nologo -MACHINE:IX86 $^ -def:$*.def \
		   $(OutPutOpt)$(EVENTLIB)
		$(LD) $(SOFLAGS) $(LDFLAGS) $^ $*.exp $(LIBS) \
		   $(OutPutOpt)$@
else
		$(LD) $(SOFLAGS) $(LDFLAGS) $^ $(OutPutOpt) $@ $(EXPLLINKLIBS)
endif
endif
	@chmod a+x $@
	@echo "done"

%.o:    %.cxx %.h
	$(CXX) $(CXXFLAGS) -c $< -o $@

clean:
	@rm -f $(OBJS) *.so G__$(PACKAGE).*

G__$(PACKAGE).cxx G__$(PACKAGE).h: $(HDRS) $(DHDR)
	@echo "Generating dictionary ..."
	rootcint -f $@ -c $(CINTFLAGS) $(ALICEINC) $^

