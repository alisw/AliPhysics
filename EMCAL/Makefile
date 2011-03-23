# $Id: Makefile 33076 2009-06-19 16:40:17Z akisiel $

include $(ROOTSYS)/test/Makefile.arch

PACKAGE = EMCAL

default-target: lib$(PACKAGE).so 

ALICEINC      = -I.

### define include dir for local case and par case

### in case of Root compilation for ANALYSIS, only load EMCALUtils
### if modification of EMCAL libraries, include STEER and RAW
ifneq ($(STEERBase_INCLUDE),)
  ALICEINC += -I../$(STEERBase_INCLUDE)
  ALICEINC += -I../$(ESD_INCLUDE)
  ALICEINC += -I../$(AOD_INCLUDE)
else
    ALICEINC += -I$(ALICE_ROOT)/include
    ALICEINC += -I$(ALICE_ROOT)/STEER
    ALICEINC += -I$(ALICE_ROOT)/RAW
    ALICEINC += -I$(ALICE_ROOT)/VZERO
endif

### in modification of EMCAL libraries, and not modification of EMCALUtils, include EMCAL from aliroot 
ifneq ($(EMCALUtils_INCLUDE),)
  ALICEINC += -I../$(EMCALUtils_INCLUDE)
else
 ALICEINC += -I$(ALICE_ROOT)/EMCAL
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

# only if no par file was loaded before
ifeq ($(ALICEINC),-I.)
  ifneq ($(ALICE_ROOT),)
    ALICEINC += -I$(ALICE_ROOT)/include
    ALICEINC += -I$(ALICE_ROOT)/STEER
    ALICEINC += -I$(ALICE_ROOT)/RAW
    ALICEINC += -I$(ALICE_ROOT)/EMCAL
  endif
endif


CXXFLAGS     += $(ALICEINC) -g
F77FLAGS     += $(ALICEINC) -g -fPIC


SRCS     = $(wildcard *.cxx */*.cxx)
HDRS     = $(SRCS:.cxx=.h)
DHDR     = $(PACKAGE)LinkDef.h

DHDR_EMCAL := $(DHDR)
HDRS_EMCAL := $(HDRS)
SRCS_EMCAL := $(SRCS) G__$(PACKAGE).cxx
FSCRS_EMCAL := $(FSRCS)
OBJS_EMCAL := $(SRCS_EMCAL:.cxx=.o) $(FSRCS:.F=.o)

PARFILE       = $(PACKAGE).par

lib$(PACKAGE).so: $(OBJS_EMCAL)
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

%.o:    %.F
	$(F77) $(F77FLAGS) -c $< -o $@

G__$(PACKAGE).cxx G__$(PACKAGE).h: $(HDRS_EMCAL) $(DHDR_EMCAL)
	@echo "Generating dictionary ..."
	rootcint -f $@ -c $(CINTFLAGS) $(ALICEINC) $^

clean:
	@rm -f $(OBJS_EMCAL) *.so G__$(PACKAGE).* $(PARFILE)

### CREATE PAR FILE

$(PARFILE): $(patsubst %,$(PACKAGE)/%,$(filter-out G__%, $(HDRS_EMCAL) $(SRCS_EMCAL) $(DHDR_EMCAL) Makefile Makefile.arch PROOF-INF))
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
