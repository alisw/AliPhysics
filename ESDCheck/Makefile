
include $(ROOTSYS)/test/Makefile.arch

default-target: libAnalysisCheck.so

ALICEINC      = -I.

### define include dir for local case and par case
ifneq ($(ANALYSIS_NEW_INCLUDE),)
  ALICEINC += -I../$(ESD_INCLUDE) -I../$(ANALYSIS_NEW_INCLUDE)
else
  ifneq ($(ALICE_ROOT),)
    ALICEINC += -I$(ALICE_ROOT)/include
  endif
endif

# for building of AnalysisCheck.par
ifneq ($(AnalysisCheck_INCLUDE),)
  ALICEINC += -I../$(AnalysisCheck_INCLUDE)
endif

CXXFLAGS     += $(ALICEINC) -g

PACKAGE = AnalysisCheck
include lib$(PACKAGE).pkg

DHDR_AnalysisCheck := $(DHDR)
HDRS_AnalysisCheck := $(HDRS)
SRCS_AnalysisCheck := $(SRCS) G__$(PACKAGE).cxx
OBJS_AnalysisCheck := $(SRCS_AnalysisCheck:.cxx=.o)

PARFILE       = $(PACKAGE).par


lib$(PACKAGE).so: $(OBJS_AnalysisCheck)
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

clean:
	@rm -f $(OBJS_AnalysisCheck) *.so G__$(PACKAGE).* $(PARFILE)

G__$(PACKAGE).cxx G__$(PACKAGE).h: $(HDRS) $(DHDR)
	@echo "Generating dictionary ..."
	rootcint -f $@ -c $(ALICEINC) $^

### CREATE PAR FILE

$(PARFILE): $(patsubst %,$(PACKAGE)/%,$(filter-out G__%, $(HDRS_AnalysisCheck) $(SRCS_AnalysisCheck) $(DHDR_AnalysisCheck) Makefile Makefile.arch lib$(PACKAGE).pkg PROOF-INF))
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
