# $Id$

PACKAGE = PWG2

include $(ROOTSYS)/test/Makefile.arch
include lib$(PACKAGE).pkg

ALICEINC      = -I. 

### define include dir for local case and par case
ifneq ($(ALICE_ROOT),)
  ALICEINC += -I$(ALICE_ROOT)/include
else
  ifneq ($(STEER_INCLUDE),)
    ALICEINC += -I../$(STEER_INCLUDE)
  endif
  ifneq ($(ESD_INCLUDE),)
    ALICEINC += -I../$(ESD_INCLUDE)
  endif
endif

CXXFLAGS     += $(ALICEINC)
SRCS         += dict.cxx
OBJS          = $(SRCS:.cxx=.o)
PARFILE       = $(PACKAGE).par

lib$(PACKAGE).so: $(OBJS)
	@echo "Linking" $@ ...
	@/bin/rm -f $@
ifeq ($(PLATFORM),macosx)
	@$(LD) -bundle -undefined $(UNDEFOPT) $(LDFLAGS) $^ -o $@
else
	@$(LD) $(SOFLAGS) $(LDFLAGS) $(OBJS) -o $@
endif
	@chmod a+x $@
	@echo "done"

%.o:    %.cxx %.h
	$(CXX) $(CXXFLAGS) -c $< -o $@

clean:
	@rm -f $(OBJS) core *~ *.so dict.* $(PARFILE)

dict.cxx dict.h: $(HDRS) $(DHDR)
	@echo "Generating dictionary ..."
	@[ -f dict ] || rm -f $@
	rootcint -f $@ -c $(ALICEINC) $^

### CREATE PAR FILE

$(PARFILE): $(patsubst %,$(PACKAGE)/%,$(filter-out dict.%, $(HDRS) $(SRCS) $(DHDR) $(PKGFILE) Makefile Makefile.arch lib$(PACKAGE).pkg PROOF-INF))
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
	@cp $^ $@

$(PACKAGE)/PROOF-INF: PROOF-INF.$(PACKAGE)
	@echo Copying $< to $@
	@[ -d $(dir $@) ] || mkdir -p $(dir $@)
	@cp -r $^ $@

$(PACKAGE)/%: %
	@echo Copying $< to $@
	@[ -d $(dir $@) ] || mkdir -p $(dir $@)
	@cp $< $@

test-%.par: %.par
	@echo "INFO: The file $< is now tested, in case of an error check in par-tmp."
	@mkdir -p par-tmp
	@cd par-tmp; tar xfz ../$<;	cd $(subst .par,,$<); PROOF-INF/BUILD.sh
	@rm -rf par-tmp
	@echo "INFO: Testing succeeded (already cleaned up)"
