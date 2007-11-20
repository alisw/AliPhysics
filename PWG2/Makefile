# $Id$

include $(ROOTSYS)/test/Makefile.arch

default-target: libPWG2.so libPWG2ebye.so libPWG2flow.so libPWG2femtoscopy.so libPWG2femtoscopyUser.so libPWG2spectra.so

ALICEINC      = -I.

### define include dir for local case and par case
ifneq ($(ALICE_ROOT),)
  ALICEINC += -I$(ALICE_ROOT)/include -I$(ALICE_ROOT)/PYTHIA6 -I$(ALICE_ROOT)/EVGEN -I$(ALICE_ROOT)/RAW -I$(ALICE_ROOT)/TPC
else
  ifneq ($(STEERBase_INCLUDE),)
    ALICEINC += -I../$(STEERBase_INCLUDE)
  endif
  ifneq ($(ESD_INCLUDE),)
    ALICEINC += -I../$(ESD_INCLUDE)
  endif
endif

# for building of PWG0dep.par
#ifneq ($(PWG0base_INCLUDE),)
#  ALICEINC += -I../$(PWG0base_INCLUDE)
#endif

CXXFLAGS     += $(ALICEINC) -g -IFEMTOSCOPY/AliFemto -IFEMTOSCOPY/AliFemto/Base -IFEMTOSCOPY/AliFemto/Infrastructure -IFEMTOSCOPY/AliFemtoUser

ALICEINC += -IFEMTOSCOPY/AliFemto -IFEMTOSCOPY/AliFemto/Base -IFEMTOSCOPY/AliFemto/Infrastructure  -IFEMTOSCOPY/AliFemtoUser

PACKAGE = PWG2
include lib$(PACKAGE).pkg

DHDR_PWG2 := $(DHDR)
HDRS_PWG2 := $(HDRS)
SRCS_PWG2 := $(SRCS) G__$(PACKAGE).cxx
OBJS_PWG2 := $(SRCS_PWG2:.cxx=.o)

PARFILE       = $(PACKAGE).par

lib$(PACKAGE).so: $(OBJS_PWG2)
	@echo "Linking" $@ ...
	@/bin/rm -f $@
ifeq ($(PLATFORM),macosx)
	@$(LD) -bundle -undefined $(UNDEFOPT) $(LDFLAGS) $^ -o $@
else
	@$(LD) $(SOFLAGS) $(LDFLAGS) $^ -o $@
endif
	@chmod a+x $@
	@echo "done"

%.o:    %.cxx %.h
	$(CXX) $(CXXFLAGS) -c $< -o $@

G__PWG2.cxx G__PWG2.h: $(HDRS_PWG2) $(DHDR_PWG2)
	@echo "Generating dictionary ..."
	rootcint -f $@ -c $(ALICEINC) $^

clean:
	@rm -f $(OBJS_PWG2) $(OBJS_PWG2ebye) $(OBJS_PWG2flow) $(OBJS_PWG2femtoscopy) *.so G__PWG2.* G__PWG2ebye.* G__PWG2flow.* G__PWG2femtoscopy.* $(PARFILE)

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

## libPWG0ebye part
include libPWG2ebye.pkg

DHDR_PWG2ebye := $(DHDR)
HDRS_PWG2ebye := $(HDRS)
SRCS_PWG2ebye := $(SRCS) G__PWG2ebye.cxx
OBJS_PWG2ebye := $(SRCS_PWG2ebye:.cxx=.o)

libPWG2ebye.so: $(OBJS_PWG2ebye)
	@echo "Linking" $@ ...
	@/bin/rm -f $@
ifeq ($(PLATFORM),macosx)
	@$(LD) -bundle -undefined $(UNDEFOPT) $(LDFLAGS) $^ -o $@
else
	@$(LD) $(SOFLAGS) $(LDFLAGS) $^ -o $@
endif
	@chmod a+x $@
	@echo "done"

G__PWG2ebye.cxx G__PWG2ebye.h: $(HDRS_PWG2ebye) $(DHDR_PWG2ebye)
	@echo "Generating dictionary ..."
	rootcint -f $@ -c $(ALICEINC) $^

## libPWG0spectra part
include libPWG2spectra.pkg

DHDR_PWG2spectra := $(DHDR)
HDRS_PWG2spectra := $(HDRS)
SRCS_PWG2spectra := $(SRCS) G__PWG2spectra.cxx
OBJS_PWG2spectra := $(SRCS_PWG2spectra:.cxx=.o)

libPWG2spectra.so: $(OBJS_PWG2spectra)
	@echo "Linking" $@ ...
	@/bin/rm -f $@
ifeq ($(PLATFORM),macosx)
	@$(LD) -bundle -undefined $(UNDEFOPT) $(LDFLAGS) $^ -o $@
else
	@$(LD) $(SOFLAGS) $(LDFLAGS) $^ -o $@
endif
	@chmod a+x $@
	@echo "done"

G__PWG2spectra.cxx G__PWG2spectra.h: $(HDRS_PWG2spectra) $(DHDR_PWG2spectra)
	@echo "Generating dictionary ..."
	rootcint -f $@ -c $(ALICEINC) $^

## libPWG2flow part
include libPWG2flow.pkg

DHDR_PWG2flow := $(DHDR)
HDRS_PWG2flow := $(HDRS)
SRCS_PWG2flow := $(SRCS) G__PWG2flow.cxx
OBJS_PWG2flow := $(SRCS_PWG2flow:.cxx=.o)

libPWG2flow.so: $(OBJS_PWG2flow)
	@echo "Linking" $@ ...
	@/bin/rm -f $@
ifeq ($(PLATFORM),macosx)
	@$(LD) -bundle -undefined $(UNDEFOPT) $(LDFLAGS) $^ -o $@
else
	@$(LD) $(SOFLAGS) $(LDFLAGS) $^ -o $@
endif
	@chmod a+x $@
	@echo "done"

G__PWG2flow.cxx G__PWG2flow.h: $(HDRS_PWG2flow) $(DHDR_PWG2flow)
	@echo "Generating dictionary ..."
	rootcint -f $@ -c $(ALICEINC) $^


## libPWG2femtoscopy part
include libPWG2femtoscopy.pkg

DHDR_PWG2femtoscopy := $(DHDR)
HDRS_PWG2femtoscopy := $(HDRS)
SRCS_PWG2femtoscopy := $(SRCS) G__PWG2femtoscopy.cxx
OBJS_PWG2femtoscopy := $(SRCS_PWG2femtoscopy:.cxx=.o)

libPWG2femtoscopy.so: $(OBJS_PWG2femtoscopy)
	@echo "Linking" $@ ...
	@/bin/rm -f $@
ifeq ($(PLATFORM),macosx)
	@$(LD) -bundle -undefined $(UNDEFOPT) $(LDFLAGS) $^ -o $@
else
	@$(LD) $(SOFLAGS) $(LDFLAGS) $^ -o $@
endif
	@chmod a+x $@
	@echo "done"

G__PWG2femtoscopy.cxx G__PWG2femtoscopy.h: $(HDRS_PWG2femtoscopy) $(DHDR_PWG2femtoscopy)
	@echo "Generating dictionary ..."
	rootcint -f $@ -c $(ALICEINC) $^

## libPWG2femtoscopyUser part
include libPWG2femtoscopyUser.pkg

DHDR_PWG2femtoscopyUser := $(DHDR)
HDRS_PWG2femtoscopyUser := $(HDRS)
SRCS_PWG2femtoscopyUser := $(SRCS) G__PWG2femtoscopyUser.cxx
OBJS_PWG2femtoscopyUser := $(SRCS_PWG2femtoscopyUser:.cxx=.o)

libPWG2femtoscopyUser.so: $(OBJS_PWG2femtoscopyUser)
	@echo "Linking" $@ ...
	@/bin/rm -f $@
ifeq ($(PLATFORM),macosx)
	@$(LD) -bundle -undefined $(UNDEFOPT) $(LDFLAGS) $^ -o $@
else
	@$(LD) $(SOFLAGS) $(LDFLAGS) $^ -o $@
endif
	@chmod a+x $@
	@echo "done"

G__PWG2femtoscopyUser.cxx G__PWG2femtoscopyUser.h: $(HDRS_PWG2femtoscopyUser) $(DHDR_PWG2femtoscopyUser)
	@echo "Generating dictionary ..."
	rootcint -f $@ -c $(ALICEINC) $^

