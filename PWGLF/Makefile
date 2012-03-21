# $Id: Makefile 50266 2011-06-24 09:37:46Z pulvir $

include $(ROOTSYS)/etc/Makefile.arch

PACKAGE = PWGLFresonancesdev

default-target: lib$(PACKAGE).so 

ALICEINC      = -I.

ifeq ($(PACKAGE),PWGLFresonances)
  ALICEINC      += -IRESONANCES
endif

ifeq ($(PACKAGE),PWGLFresonancesdev)
  ALICEINC      += -IRESONANCES_DEV
endif

### define include dir for local case and par case
ifneq ($(STEERBase_INCLUDE),)
   ALICEINC += -I../$(STEERBase_INCLUDE)
   ALICEINC += -I../$(ESD_INCLUDE)
   ALICEINC += -I../$(AOD_INCLUDE)
   ALICEINC += -I../$(ANALYSIS_INCLUDE)
   ALICEINC += -I../$(ANALYSISalice_INCLUDE)
else
   ifneq ($(ALICE_ROOT),)
      ALICEINC += -I$(ALICE_ROOT)/include
   endif
endif

ifneq ($(CORRFW_INCLUDE),)
   ALICEINC += -I../$(CORRFW_INCLUDE)
endif

ifneq ($(PWGLFresonances_INCLUDE),)
   ALICEINC += -I../$(PWGLFresonances_INCLUDE)
else
   ifneq ($(ALICE_ROOT),)
      ALICEINC += -I$(ALICE_ROOT)/PWGLF -I$(ALICE_ROOT)/PWGLF/RESONANCES
   endif
endif

ifneq ($(PWGLFresonancesdev_INCLUDE),)
   ALICEINC += -I../$(PWGLFresonancesdev_INCLUDE)
else
   ifneq ($(ALICE_ROOT),)
      ALICEINC += -I$(ALICE_ROOT)/PWGLF -I$(ALICE_ROOT)/PWGLF/RESONANCES_DEV
   endif
endif

# only if no par file was loaded before
ifeq ($(ALICEINC),-I.)
  ifneq ($(ALICE_ROOT),)
    ALICEINC += -I$(ALICE_ROOT)/include
  endif
endif

ifneq ($(EINCLUDE),)
  ALICEINC += $(EINCLUDE:%=-I$(ALICE_ROOT)/%)
endif

CXXFLAGS     += $(ALICEINC) -g
F77FLAGS     += $(ALICEINC) -g -fPIC

DHDR = 
HDRS = 
SRCS = 
FSRCS = 

DHDR_PWG2 := $(DHDR)
HDRS_PWG2 := $(HDRS)
SRCS_PWG2 := $(SRCS) G__$(PACKAGE).cxx
FSCRS_PWG2 := $(FSRCS)
OBJS_PWG2 := $(SRCS_PWG2:.cxx=.o) $(FSRCS:.F=.o)

PARFILE       = $(PACKAGE).par

lib$(PACKAGE).so: $(OBJS_PWG2)
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

G__$(PACKAGE).h:G__$(PACKAGE).cxx
	@if test -f $@ ; then : ; else rm -f $< ; $(MAKE) $@ ; fi

G__$(PACKAGE).cxx: $(HDRS_PWG2) $(DHDR_PWG2)
	@echo "Generating dictionary ..."
	rootcint -f $@ -c $(CINTFLAGS) $(ALICEINC) $^

clean:
	@rm -f $(OBJS_PWG2) *.so G__$(PACKAGE).* $(PARFILE)

### CREATE PAR FILE

$(PARFILE): $(patsubst %,$(PACKAGE)/%,$(filter-out G__%, $(HDRS_PWG2) $(SRCS_PWG2) $(DHDR_PWG2) Makefile Makefile.arch PROOF-INF))
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
