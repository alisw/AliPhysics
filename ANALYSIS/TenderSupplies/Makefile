
include $(ROOTSYS)/test/Makefile.arch

PACKAGE = TenderSupplies

default-target: lib$(PACKAGE).so

ALICEINC      = -I.

# add include paths from other par files
ifneq ($(ESD_INCLUDE),)
	ALICEINC += -I../$(ESD_INCLUDE)
endif
	
ifneq ($(AOD_INCLUDE),)
	ALICEINC += -I../$(AOD_INCLUDE)
endif

ifneq ($(STEERBase_INCLUDE),)
	ALICEINC += -I../$(STEERBase_INCLUDE)
endif
			
ifneq ($(ANALYSIS_INCLUDE),)
	ALICEINC += -I../$(ANALYSIS_INCLUDE)
endif
				
ifneq ($(ANALYSISalice_INCLUDE),)
	ALICEINC += -I../$(ANALYSISalice_INCLUDE)
endif

		
# only if no par file was loaded before
ifeq ($(ALICEINC),-I.)
	ifneq ($(ALICE_ROOT),)
		ALICEINC += -I./dielectron -I$(ALICE_ROOT)/include -I$(ALICE_ROOT)/STEER -I$(ALICE_ROOT)/ANALYSIS -I$(ALICE_ROOT)/ANALYSIS/Tender -I$(ALICE_ROOT)/TPC -I$(ALICE_ROOT)/TOF
	endif
endif

#include HFE in case of the train use the train checkout
ifneq ($(TRAIN_ROOT),)
# 	ALICEINC += -I$(TRAIN_ROOT)/util/hfe/hfe
	CXXFLAGS     += $(ALICEINC) -g
else
# 	ALICEINC += -I$(ALICE_ROOT)/PWG3/hfe
	CXXFLAGS     += -W -Wall -Weffc++ -Woverloaded-virtual -fPIC -pipe -fmessage-length=0 -Wno-long-long -pedantic-errors -ansi -Dlinux $(ALICEINC) -g
endif
				
		
include lib$(PACKAGE).pkg
						
DHDR_TenderSupplies := $(DHDR)
HDRS_TenderSupplies := $(HDRS)
SRCS_TenderSupplies := $(SRCS) G__$(PACKAGE).cxx
OBJS_TenderSupplies := $(SRCS_TenderSupplies:.cxx=.o)
						
PARFILE       = $(PACKAGE).par
						
						
lib$(PACKAGE).so: $(OBJS_TenderSupplies)
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
	$(CXX) $(CXXFLAGS) $(PACKCXXFLAGS) -c $< -o $@

clean:
	@rm -f $(OBJS_TenderSupplies) *.so G__$(PACKAGE).* $(PARFILE)

G__$(PACKAGE).cxx G__$(PACKAGE).h: $(HDRS) $(DHDR)
	@echo "Generating dictionary ..."
	rootcint -f $@ -c $(CINTFLAGS) $(ALICEINC) $^

### CREATE PAR FILE

$(PARFILE): $(patsubst %,$(PACKAGE)/%,$(filter-out G__%, $(HDRS_TenderSupplies) $(SRCS_TenderSupplies) $(DHDR_TenderSupplies) Makefile Makefile.arch lib$(PACKAGE).pkg PROOF-INF))
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
	
