# $Id: Makefile 53653 2011-12-16 10:05:45Z cholm $

include Makefile.arch

PACKAGE = PWGLFresonances

default-target: lib$(PACKAGE).so 

ALICEINC      = -I. -I$(ALICE_ROOT)/TOF -I$(ALICE_ROOT)/ANALYSIS/EventMixing

### define include dir for local case and par case
ifneq ($(ESD_INCLUDE),)
   ALICEINC += -I../$(ESD_INCLUDE)
   ALICEINC += -I../$(AOD_INCLUDE)
   ALICEINC += -I../$(STEERBase_INCLUDE)
   ALICEINC += -I../$(OADB_INCLUDE)
   ALICEINC += -I../$(ANALYSIS_INCLUDE)
   ALICEINC += -I../$(ANALYSISalice_INCLUDE)
 else
    ifneq ($(ALICE_ROOT),)
       ALICEINC += -I$(ALICE_ROOT)/include -I$(ALICE_ROOT)/PWGLF
     endif
endif

ifneq ($(PWGLFresonances_INCLUDE),)
  ALICEINC += -I../$(PWGLFresonances_INCLUDE)
 else
    ifneq ($(ALICE_ROOT),)
       ALICEINC += -I$(ALICE_ROOT)/PWGLF/RESONANCES
     endif
endif

ifneq ($(PWGLFSTRANGENESS_INCLUDE),)
  ALICEINC += -I../$(PWGLFSTRANGENESS_INCLUDE)
 else
    ifneq ($(ALICE_ROOT),)
       ALICEINC += -I$(ALICE_ROOT)/PWGLF/STRANGENESS
     endif
endif

ifneq ($(PWGLFtotEt_INCLUDE),)
  ALICEINC += -I../$(PWGLFtotEt_INCLUDE)
 else
    ifneq ($(ALICE_ROOT),)
       ALICEINC += -I$(ALICE_ROOT)/PWGLF/totEt
     endif
endif

ifneq ($(PWGLFforward_INCLUDE),)
  ALICEINC += -I../$(PWGLFforward_INCLUDE)
 else
    ifneq ($(ALICE_ROOT),)
       ALICEINC += -I$(ALICE_ROOT)/PWGLF/FORWARD/analysis
     endif
endif

ifneq ($(PWGLFforward2_INCLUDE),)
  ALICEINC += -I../$(PWGLFforward2_INCLUDE)
 else
    ifneq ($(ALICE_ROOT),)
       ALICEINC += -I$(ALICE_ROOT)/PWGLF/FORWARD/analysis2
     endif
endif

ifneq ($(OADB_INCLUDE), )
  ALICEINC += -I../$(OADB_INCLUDE)
else
  ifneq ($(ALICE_ROOT),)
    ALICEINC += -I$(ALICE_ROOT)/OADB
  endif
endif

# only if no par file was loaded before
#ifeq ($(ALICEINC),-I.)
  ifneq ($(ALICE_ROOT),)
    ALICEINC += -I$(ALICE_ROOT)/include
  endif
#endif

ifneq ($(EINCLUDE),)
  ALICEINC += $(EINCLUDE:%=-I$(ALICE_ROOT)/%)
endif

CXXFLAGS     += $(ALICEINC) -g
F77FLAGS     += $(ALICEINC) -g -fPIC

DHDR =
HDRS =
SRCS =
FSRCS =

DHDR_PWGLF := $(DHDR)
HDRS_PWGLF := $(HDRS)
SRCS_PWGLF := $(SRCS) G__$(PACKAGE).cxx
FSCRS_PWGLF := $(FSRCS)
OBJS_PWGLF := $(SRCS_PWGLF:.cxx=.o) $(FSRCS:.F=.o)

PARFILE       = $(PACKAGE).par

lib$(PACKAGE).so: $(OBJS_PWGLF)
	@echo "Linking" $@ ...
	@/bin/rm -f $@
	$(LD) $(SOFLAGS) $(LDFLAGS) $^ $(OutPutOpt) $@ $(EXPLLINKLIBS)
	@chmod a+x $@
	@echo "done"

%.o:    %.cxx %.h
	$(CXX) $(CXXFLAGS) -c $< -o $@

%.o:    %.F
	$(F77) $(F77FLAGS) -c $< -o $@

G__$(PACKAGE).h:G__$(PACKAGE).cxx
	@if test -f $@ ; then : ; else rm -f $< ; $(MAKE) $@ ; fi

G__$(PACKAGE).cxx: $(HDRS_PWGLF) $(DHDR_PWGLF)
	@echo "Generating dictionary ..."
	rootcint -f $@ -c $(CINTFLAGS) $(ALICEINC) $^

clean:
	@rm -f $(OBJS_PWGLF) *.so G__$(PACKAGE).* $(PARFILE)

### CREATE PAR FILE

$(PARFILE): $(patsubst %,$(PACKAGE)/%,$(filter-out G__%, $(HDRS_PWGLF) $(SRCS_PWGLF) $(DHDR_PWGLF) Makefile Makefile.arch lib$(PACKAGE).pkg PROOF-INF))
	@echo "Creating archive" $@ ...
	@tar cfzh $@ $(PACKAGE)
	@rm -rf $(PACKAGE)
	@echo "done"

$(PACKAGE)/Makefile: Makefile #.$(PACKAGE)
	@echo Copying $< to $@ with transformations
	@[ -d $(dir $@) ] || mkdir -p $(dir $@)
	@sed 's/include \$$(ROOTSYS)\/test\/Makefile.arch/include Makefile.arch/' < $^ > $@

$(PACKAGE)/Makefile.arch: $(ROOTSYS)/test/Makefile.arch
	@echo Copying $< to $@
	@[ -d $(dir $@) ] || mkdir -p $(dir $@)
	@cp -a $^ $@

$(PACKAGE)/PROOF-INF: PROOF-INF.$(PACKAGE)
	@echo Copying $< to $@
	@[ -d $(dir $@) ] || mkdir -p $(dir $@)
	@cp -a -r $^ $@

$(PACKAGE)/%: %
	@echo Copying $< to $@
	@[ -d $(dir $@) ] || mkdir -p $(dir $@)
	@cp -a $< $@

test-%.par: %.par
	@echo "INFO: The file $< is now tested, in case of an error check in par-tmp."
	@mkdir -p par-tmp
	@cd par-tmp; tar xfz ../$<;	cd $(subst .par,,$<); PROOF-INF/BUILD.sh
	@rm -rf par-tmp
	@echo "INFO: Testing succeeded (already cleaned up)"
