# $Id$

include $(ROOTSYS)/test/Makefile.arch

default-target: libPWG0base.so libPWG0dep.so libPWG0selectors.so

ALICEINC      = -I.

### define include dir for local case and par case
ifneq ($(ALICE_ROOT),)
  ALICEINC += -I$(ALICE_ROOT)/include -I$(ALICE_ROOT)/PYTHIA6 -I$(ALICE_ROOT)/EVGEN -I$(ALICE_ROOT)/RAW -I$(ALICE_ROOT)/TPC
else
  ifneq ($(ESD_INCLUDE),)
    ALICEINC += -I../$(ESD_INCLUDE)
  endif
endif

# for building of PWG0dep.par
ifneq ($(PWG0base_INCLUDE),)
  ALICEINC += -I../$(PWG0base_INCLUDE)
endif

CXXFLAGS     += $(ALICEINC) -g

PACKAGE = PWG0base
include lib$(PACKAGE).pkg

DHDR_PWG0base := $(DHDR)
HDRS_PWG0base := $(HDRS)
SRCS_PWG0base := $(SRCS) G__$(PACKAGE).cxx
OBJS_PWG0base := $(SRCS_PWG0base:.cxx=.o)

PARFILE       = $(PACKAGE).par

lib$(PACKAGE).so: $(OBJS_PWG0base)
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

G__PWG0base.cxx G__PWG0base.h: $(HDRS_PWG0base) $(DHDR_PWG0base)
	@echo "Generating dictionary ..."
	rootcint -f $@ -c $(ALICEINC) $^

clean:
	@rm -f $(OBJS_PWG0base) $(OBJS_PWG0dep) $(OBJS_PWG0selectors) *.so G__PWG0base.* G__PWG0dep.* G__PWG0selectors.* $(PARFILE)

### CREATE PAR FILE

$(PARFILE): $(patsubst %,$(PACKAGE)/%,$(filter-out G__%, $(HDRS_PWG0base) $(SRCS_PWG0base) $(DHDR_PWG0base) Makefile Makefile.arch lib$(PACKAGE).pkg PROOF-INF))
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

## libPWG0dep part
include libPWG0dep.pkg

DHDR_PWG0dep := $(DHDR)
HDRS_PWG0dep := $(HDRS)
SRCS_PWG0dep := $(SRCS) G__PWG0dep.cxx
OBJS_PWG0dep := $(SRCS_PWG0dep:.cxx=.o)

libPWG0dep.so: $(OBJS_PWG0dep)
	@echo "Linking" $@ ...
	@/bin/rm -f $@
ifeq ($(PLATFORM),macosx)
	@$(LD) -bundle -undefined $(UNDEFOPT) $(LDFLAGS) $^ -o $@
else
	@$(LD) $(SOFLAGS) $(LDFLAGS) $^ -o $@
endif
	@chmod a+x $@
	@echo "done"

G__PWG0dep.cxx G__PWG0dep.h: $(HDRS_PWG0dep) $(DHDR_PWG0dep)
	@echo "Generating dictionary ..."
	rootcint -f $@ -c $(ALICEINC) $^

## libPWG0selectors part
include libPWG0selectors.pkg

DHDR_PWG0selectors := $(DHDR)
HDRS_PWG0selectors := $(HDRS)
SRCS_PWG0selectors := $(SRCS) G__PWG0selectors.cxx
OBJS_PWG0selectors := $(SRCS_PWG0selectors:.cxx=.o)

libPWG0selectors.so: $(OBJS_PWG0selectors)
	@echo "Linking" $@ ...
	@/bin/rm -f $@
ifeq ($(PLATFORM),macosx)
	@$(LD) -bundle -undefined $(UNDEFOPT) $(LDFLAGS) $^ -o $@
else
	@$(LD) $(SOFLAGS) $(LDFLAGS) $^ -o $@
endif
	@chmod a+x $@
	@echo "done"

G__PWG0selectors.cxx G__PWG0selectors.h: $(HDRS_PWG0selectors) $(DHDR_PWG0selectors)
	@echo "Generating dictionary ..."
	rootcint -f $@ -c $(ALICEINC) $^

