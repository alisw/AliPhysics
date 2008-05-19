
include Makefile.arch

default-target: libPWG4PartCorr.so

ALICEINC      = -I.

### define include dir for local case and par case
ifneq ($(ANALYSIS_INCLUDE),)
  ALICEINC += -I../$(ESD_INCLUDE) -I../$(AOD_INCLUDE) -I../$(ANALYSIS_INCLUDE) -I../$(STEERBase_INCLUDE) -I../$(ANALYSISalice_INCLUDE)
else
  ifneq ($(ALICE_ROOT),)
    ALICEINC += -I$(ALICE_ROOT)/include
  endif
endif

# for building of PWG4PartCorr.par
ifneq ($(GAMMA_INCLUDE),)
  ALICEINC += -I../$(GAMMA_INCLUDE)
endif

CXXFLAGS     += $(ALICEINC) -g

PACKAGE = PWG4PartCorr
include lib$(PACKAGE).pkg

DHDR_PWG4PartCorr := $(DHDR)
HDRS_PWG4PartCorr := $(HDRS)
SRCS_PWG4PartCorr := $(SRCS) G__$(PACKAGE).cxx
OBJS_PWG4PartCorr := $(SRCS_PWG4PartCorr:.cxx=.o)

PARFILE       = $(PACKAGE).par


lib$(PACKAGE).so: $(OBJS_PWG4PartCorr)
	@echo "Linking" $@ ...
	@/bin/rm -f $@
ifeq ($(ARCH),macosx)
	@$(LD) -bundle -undefined $(UNDEFOPT) $(LDFLAGS) $^ -o $@
else
	@$(LD) $(SOFLAGS) $(LDFLAGS) $^ -o $@
endif
	@chmod a+x $@
	@echo "done"

%.o:    %.cxx %.h
	$(CXX) $(CXXFLAGS) -c $< -o $@

clean:
	@rm -f $(OBJS_PWG4PartCorr) *.so G__$(PACKAGE).* $(PARFILE)

G__$(PACKAGE).cxx G__$(PACKAGE).h: $(HDRS) $(DHDR)
	@echo "Generating dictionary ..."
	rootcint -f $@ -c $(ALICEINC) $^

### CREATE PAR FILE

$(PARFILE): $(patsubst %,$(PACKAGE)/%,$(filter-out G__%, $(HDRS_PWG4PartCorr) $(SRCS_PWG4PartCorr) $(DHDR_PWG4PartCorr) Makefile Makefile.arch lib$(PACKAGE).pkg PROOF-INF))
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
