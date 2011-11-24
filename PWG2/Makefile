# $Id$

include $(ROOTSYS)/test/Makefile.arch

PACKAGE = PWG2resonances

default-target: lib$(PACKAGE).so 

ALICEINC      = -I. -I$(ALICE_ROOT)/TOF -I$(ALICE_ROOT)/ANALYSIS/EventMixing

### define include dir for local case and par case
ifneq ($(ESD_INCLUDE),)
   ALICEINC += -I../$(ESD_INCLUDE)
   ALICEINC += -I../$(AOD_INCLUDE)
   ALICEINC += -I../$(STEERBase_INCLUDE)
   ALICEINC += -I../$(ANALYSIS_INCLUDE)
   ALICEINC += -I../$(ANALYSISalice_INCLUDE)
 else
    ifneq ($(ALICE_ROOT),)
       ALICEINC += -I$(ALICE_ROOT)/include -I$(ALICE_ROOT)/PWG2
     endif
endif

ifneq ($(CORRFW_INCLUDE),)
  ALICEINC += -I../$(CORRFW_INCLUDE)
endif

ifneq ($(PWG2AOD_INCLUDE),)
  ALICEINC += -I../$(PWG2AOD_INCLUDE)
endif

ifneq ($(PWG2femtoscopy_INCLUDE),)
  ALICEINC += -I../$(PWG2femtoscopy_INCLUDE)
 else
    ifneq ($(ALICE_ROOT),)
       ALICEINC += -I$(ALICE_ROOT)/PWG2/AOD
     endif
endif

ifneq ($(PWG2femtoscopyUser_INCLUDE),)
  ALICEINC += -I../$(PWG2femtoscopyUser_INCLUDE)
 else
    ifneq ($(ALICE_ROOT),)
       ALICEINC += -I$(ALICE_ROOT)/PWG2/AOD -I$(ALICE_ROOT)/PWG2/FEMTOSCOPY/AliFemto
     endif
endif

ifneq ($(PWG2flowCommon_INCLUDE),)
  ALICEINC += -I../$(PWG2flowCommon_INCLUDE)
endif

ifneq ($(PWG2flowTasks_INCLUDE),)
  ALICEINC += -I../$(PWG2flowTasks_INCLUDE)
 else
    ifneq ($(ALICE_ROOT),)
       ALICEINC += -I$(ALICE_ROOT)/PWG2/FLOW/AliFlowCommon
     endif
endif

ifneq ($(PWG2resonances_INCLUDE),)
  ALICEINC += -I../$(PWG2resonances_INCLUDE)
 else
    ifneq ($(ALICE_ROOT),)
       ALICEINC += -I$(ALICE_ROOT)/PWG2/RESONANCES
     endif
endif

ifneq ($(PWG2unicor_INCLUDE),)
  ALICEINC += -I../$(PWG2unicor_INCLUDE)
endif

ifneq ($(PWG2evchar_INCLUDE),)
  ALICEINC += -I../$(PWG2evchar_INCLUDE)
endif

ifneq ($(PWG2kink_INCLUDE),)
  ALICEINC += -I../$(PWG2kink_INCLUDE)
endif

ifneq ($(PWG2forward_INCLUDE),)
  ALICEINC += -I../$(PWG2forward_INCLUDE)
 else
    ifneq ($(ALICE_ROOT),)
       ALICEINC += -I$(ALICE_ROOT)/FMD
     endif
endif

# only if no par file was loaded before
ifeq ($(ALICEINC),-I.)
  ifneq ($(ALICE_ROOT),)
    ALICEINC += -I$(ALICE_ROOT)/include
  endif
endif

# for building of PWG0dep.par
#ifneq ($(PWG0base_INCLUDE),)
#  ALICEINC += -I../$(PWG0base_INCLUDE)
#endif

#ALICEINC += -IFEMTOSCOPY/AliFemto -IANALYSIS -IFEMTOSCOPY/AliFemtoUser

# include lib$(PACKAGE).pkg

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

$(PARFILE): $(patsubst %,$(PACKAGE)/%,$(filter-out G__%, $(HDRS_PWG2) $(SRCS_PWG2) $(DHDR_PWG2) Makefile Makefile.arch lib$(PACKAGE).pkg PROOF-INF))
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
