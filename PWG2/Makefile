# $Id$

include $(ROOTSYS)/test/Makefile.arch

default-target: libPWG2.so libPWG2ebye.so libPWG0flow.so 

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
#ifneq ($(PWG0base_INCLUDE),)
#  ALICEINC += -I../$(PWG0base_INCLUDE)
#endif

CXXFLAGS     += $(ALICEINC) -g

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
	@rm -f $(OBJS_PWG2) $(OBJS_PWG2ebye) $(OBJS_PWG2flow) *.so G__PWG2.* G__PWG2ebye.* G__PWG2flow.* $(PARFILE)

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

