# $Id$

include $(ROOTSYS)/test/Makefile.arch

default-target: libPWG3base.so libPWG3muon.so libPWG3muondep.so libPWG3vertexingHF.so libPWG3hfe.so libPWG3dielectron.so

ALICEINC      = -I.

### define include dir for local case and par case
ifneq ($(ESD_INCLUDE),)
   ALICEINC += -I../$(ESD_INCLUDE)
	ifneq ($(STEERBase_INCLUDE),)
		ALICEINC += -I../$(STEERBase_INCLUDE)
	endif				
	ifneq ($(AOD_INCLUDE),)
		ALICEINC += -I../$(AOD_INCLUDE)
	endif
	ifneq ($(ANALYSIS_INCLUDE),)
		ALICEINC += -I../$(ANALYSIS_INCLUDE)
	endif
	ifneq ($(ANALYSISalice_INCLUDE),)
		ALICEINC += -I../$(ANALYSISalice_INCLUDE)
	endif
	ifneq ($(CORRFW_INCLUDE),)
		ALICEINC += -I../$(CORRFW_INCLUDE)
	endif
else
  ifneq ($(ALICE_ROOT),)
    ALICEINC += -I$(ALICE_ROOT)/include -I$(ALICE_ROOT)/MUON -I$(ALICE_ROOT)/MUON/mapping
else
  endif
endif

CXXFLAGS     += $(ALICEINC) -g

PACKAGE = PWG3
include lib$(PACKAGE).pkg

DHDR_PWG3 := $(DHDR)
HDRS_PWG3 := $(HDRS)
SRCS_PWG3 := $(SRCS) G__$(PACKAGE).cxx
OBJS_PWG3 := $(SRCS_PWG3:.cxx=.o)

PARFILE       = $(PACKAGE).par

lib$(PACKAGE).so: $(OBJS_PWG3)
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

G__PWG3.cxx G__PWG3.h: $(HDRS_PWG3) $(DHDR_PWG3)
	@echo "Generating dictionary ..."
	rootcint -f $@ -c $(ALICEINC) $^

clean:
	@rm -f $(OBJS_PWG3)  *.so G__PWG3.*  $(PARFILE)

### CREATE PAR FILE

$(PARFILE): $(patsubst %,$(PACKAGE)/%,$(filter-out G__%, $(HDRS_PWG3) $(SRCS_PWG3) $(DHDR_PWG3) Makefile Makefile.arch lib$(PACKAGE).pkg PROOF-INF))
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

## libPWG3base part
##include libPWG3base.pkg

DHDR_PWG3base := $(DHDR)
HDRS_PWG3base := $(HDRS)
SRCS_PWG3base := $(SRCS) G__PWG3base.cxx
OBJS_PWG3base := $(SRCS_PWG3base:.cxx=.o)

libPWG3base.so: $(OBJS_PWG3base)
	@echo "Linking" $@ ...
	@/bin/rm -f $@
ifeq ($(ARCH),macosx)
	@$(LD) -bundle -undefined $(UNDEFOPT) $(LDFLAGS) $^ -o $@
else
	@$(LD) $(SOFLAGS) $(LDFLAGS) $^ -o $@
endif
	@chmod a+x $@
	@echo "done"

G__PWG3base.cxx G__PWG3base.h: $(HDRS_PWG3base) $(DHDR_PWG3base)
	@echo "Generating dictionary ..."
	rootcint -f $@ -c $(ALICEINC) $^


## libPWG3muon part
##include libPWG3muon.pkg

DHDR_PWG3muon := $(DHDR)
HDRS_PWG3muon := $(HDRS)
SRCS_PWG3muon := $(SRCS) G__PWG3muon.cxx
OBJS_PWG3muon := $(SRCS_PWG3muon:.cxx=.o)

libPWG3muon.so: $(OBJS_PWG3muon)
	@echo "Linking" $@ ...
	@/bin/rm -f $@
ifeq ($(ARCH),macosx)
	@$(LD) -bundle -undefined $(UNDEFOPT) $(LDFLAGS) $^ -o $@
else
	@$(LD) $(SOFLAGS) $(LDFLAGS) $^ -o $@
endif
	@chmod a+x $@
	@echo "done"

G__PWG3muon.cxx G__PWG3muon.h: $(HDRS_PWG3muon) $(DHDR_PWG3muon)
	@echo "Generating dictionary ..."
	rootcint -f $@ -c $(ALICEINC) $^

## libPWG3muondep part
##include libPWG3muondep.pkg

DHDR_PWG3muondep := $(DHDR)
HDRS_PWG3muondep := $(HDRS)
SRCS_PWG3muondep := $(SRCS) G__PWG3muondep.cxx
OBJS_PWG3muondep := $(SRCS_PWG3muondep:.cxx=.o)

libPWG3muondep.so: $(OBJS_PWG3muondep)
	@echo "Linking" $@ ...
	@/bin/rm -f $@
ifeq ($(ARCH),macosx)
	@$(LD) -bundle -undefined $(UNDEFOPT) $(LDFLAGS) $^ -o $@
else
	@$(LD) $(SOFLAGS) $(LDFLAGS) $^ -o $@
endif
	@chmod a+x $@
	@echo "done"

G__PWG3muondep.cxx G__PWG3muondep.h: $(HDRS_PWG3muondep) $(DHDR_PWG3muondep)
	@echo "Generating dictionary ..."
	rootcint -f $@ -c $(ALICEINC) $^

## libPWG3vertexingHF part
##include libPWG3vertexingHF.pkg

DHDR_PWG3vertexingHF := $(DHDR)
HDRS_PWG3vertexingHF := $(HDRS)
SRCS_PWG3vertexingHF := $(SRCS) G__PWG3vertexingHF.cxx
OBJS_PWG3vertexingHF := $(SRCS_PWG3vertexingHF:.cxx=.o)

libPWG3vertexingHF.so: $(OBJS_PWG3vertexingHF)
	@echo "Linking" $@ ...
	@/bin/rm -f $@
ifeq ($(ARCH),macosx)
	@$(LD) -bundle -undefined $(UNDEFOPT) $(LDFLAGS) $^ -o $@
else
	@$(LD) $(SOFLAGS) $(LDFLAGS) $^ -o $@
endif
	@chmod a+x $@
	@echo "done"

G__PWG3vertexingHF.cxx G__PWG3vertexingHF.h: $(HDRS_PWG3vertexingHF) $(DHDR_PWG3vertexingHF)
	@echo "Generating dictionary ..."
	rootcint -f $@ -c $(ALICEINC) $^


## libPWG3hfe part
##include libPWG3hfe.pkg

DHDR_PWG3hfe := $(DHDR)
HDRS_PWG3hfe := $(HDRS)
SRCS_PWG3hfe := $(SRCS) G__PWG3hfe.cxx
OBJS_PWG3hfe := $(SRCS_PWG3hfe:.cxx=.o)

libPWG3hfe.so: $(OBJS_PWG3hfe)
	@echo "Linking" $@ ...
	@/bin/rm -f $@
ifeq ($(ARCH),macosx)
	@$(LD) -bundle -undefined $(UNDEFOPT) $(LDFLAGS) $^ -o $@
else
	@$(LD) $(SOFLAGS) $(LDFLAGS) $^ -o $@
endif
	@chmod a+x $@
	@echo "done"

G__PWG3hfe.cxx G__PWG3hfe.h: $(HDRS_PWG3hfe) $(DHDR_PWG3hfe)
	@echo "Generating dictionary ..."
	rootcint -f $@ -c $(ALICEINC) $^

## libPWG3dielectron part
##include libPWG3dielectron.pkg
	
DHDR_PWG3dielectron := $(DHDR)
HDRS_PWG3dielectron := $(HDRS)
SRCS_PWG3dielectron := $(SRCS) G__PWG3dielectron.cxx
OBJS_PWG3dielectron := $(SRCS_PWG3dielectron:.cxx=.o)

libPWG3dielectron.so: $(OBJS_PWG3dielectron)
	@echo "Linking" $@ ...
	@/bin/rm -f $@
ifeq ($(ARCH),macosx)
	@$(LD) -bundle -undefined $(UNDEFOPT) $(LDFLAGS) $^ -o $@
else
	@$(LD) $(SOFLAGS) $(LDFLAGS) $^ -o $@
endif
	@chmod a+x $@
	@echo "done"
			
G__PWG3dielectron.cxx G__PWG3dielectron.h: $(HDRS_PWG3dielectron) $(DHDR_PWG3dielectron)
	@echo "Generating dictionary ..."
	rootcint -f $@ -c $(ALICEINC) $^
