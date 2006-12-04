# $Id$

include $(ROOTSYS)/test/Makefile.arch

default-target: libPWG2ebye.so

ALICEINC      = -I.

### define include dir for local case and par case
ifneq ($(ALICE_ROOT),)
  ALICEINC += -I$(ALICE_ROOT)/include -I$(ALICE_ROOT)/PYTHIA6 -I$(ALICE_ROOT)/EVGEN
else
  ifneq ($(ESD_INCLUDE),)
    ALICEINC += -I../$(ESD_INCLUDE)
  endif
endif

CXXFLAGS     += $(ALICEINC) -g

PACKAGE = PWG2ebye
include lib$(PACKAGE).pkg

DHDR_PWG2ebye := $(DHDR)
HDRS_PWG2ebye := $(HDRS)
SRCS_PWG2ebye := $(SRCS) G__$(PACKAGE).cxx
OBJS_PWG2ebye := $(SRCS_PWG2ebye:.cxx=.o)

PARFILE       = $(PACKAGE).par

lib$(PACKAGE).so: $(OBJS_PWG2ebye)
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

G__PWG2ebye.cxx G__PWG2ebye.h: $(HDRS_PWG2ebye) $(DHDR_PWG2ebye)
	@echo "Generating dictionary ..."
	rootcint -f $@ -c $(ALICEINC) $^

clean:
	@rm -f $(OBJS_PWG2ebye) *.so G__PWG2ebye.* $(PARFILE)

### CREATE PAR FILE

$(PARFILE): $(patsubst %,$(PACKAGE)/%,$(filter-out G__%, $(HDRS_PWG2ebye) $(SRCS_PWG2ebye) $(DHDR_PWG2ebye) Makefile Makefile.arch lib$(PACKAGE).pkg PROOF-INF))
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

