############################### Main Makefile #################################

# Include machine specific definitions

include $(ALICE_ROOT)/conf/GeneralDef
include $(ALICE_ROOT)/conf/MachineDef.$(ALICE_TARGET)

MAKEFLAGS += -s

##### MACROS #####

PACKAGE = Main

DOTS = " ................................................................................"

PRETTY =  $(AWK) '{print $$0 substr($(DOTS),1,79-length($$0))}'

##### Module libraries #####

ALIROOT_DIRS		= STEER TGeant3 TRD PHOS TPC ZDC MUON PMD FMD TOF ITS \
			  CASTOR RICH START STRUCT EVGEN RALICE ALIFAST THijing

##### TARGETS #####
 
default:      lib bin alilibs aliroot

lib bin:
	@mkdir $@

alilibs:  lib
	echo MAKEFLAGS = $(MAKEFLAGS)
	for i in $(ALIROOT_DIRS) ; do \
	   echo "Making headers in $$i" | $(PRETTY); \
	   ${MAKE} -C $$i headers ; \
        done
	@for i in $(ALIROOT_DIRS) ; do \
	   echo "Making dependencies in $$i" | $(PRETTY); \
	   ${MAKE} -C $$i depend ; \
	done
	@for i in $(ALIROOT_DIRS) ; do \
	   echo "Making in $$i" | $(PRETTY); \
	   ${MAKE} -C $$i ; \
	done

aliroot geant321 minicern pdf pythia hijing: FORCE
	@DIR=`echo $@ | $(AWK) '{print toupper($$0)}'` ; \
	echo "Making dependencies in $$DIR" | $(PRETTY); \
	${MAKE} -C $$DIR depend;\
	echo "Making in $$DIR" | $(PRETTY); \
	${MAKE} -C $$DIR

TGeant4 AliGeant4: FORCE
	@DIR=$@; \
	echo "Making dependencies in $$DIR" | $(PRETTY); \
	${MAKE} -C $$DIR depend;\
	echo "Making in $$DIR" | $(PRETTY); \
	${MAKE} -C $$DIR

cernlibs: geant321 pythia minicern pdf hijing

geant4: TGeant4 AliGeant4

all:	cernlibs default

FORCE:

############################### General Macros ################################

# include $(ALICE_ROOT)/conf/GeneralMacros

############################### Specific Macros ###############################

STRUCT_DIRS	= html conf macros data share include Euclid picts \
                  doc etc Makefile .rootrc

LIBRARY_DIRS	= MINICERN GEANT321 PYTHIA PDF

dist: AliRoot$(VERSION).tar.gz

AliRoot$(VERSION).tar.gz: $(STRUCT_DIRS) $(ALIROOT_DIRS) ALIROOT

distall: AliOffline$(VERSION).tar.gz

AliOffline$(VERSION).tar.gz: $(STRUCT_DIRS) $(ALIROOT_DIRS) $(LIBRARY_DIRS) ALIROOT

distlib: AliLibs$(VERSION).tar.gz

AliLibs$(VERSION).tar.gz: $(LIBRARY_DIRS)

AliRoot$(VERSION).tar.gz AliLibs$(VERSION).tar.gz AliOffline$(VERSION).tar.gz:
		@rm -f $(ALICE)/$@ 
		@rm -f `find . -name '*~' -print` \
                       `find . -name '*.bak' -print` \
                       `find . -name '.*~' -print` \
		       `find . -name '*\#*' -print` 
		@rm -f /tmp/saves
		@ls -1d $^ | sed -e "s/^/$(ALICE_LEVEL)\//" > /tmp/saves
	 	@cd $(ALICE) ; \
                gtar cvfz $@ --exclude '*.o' --exclude '*Cint.*' \
                --exclude 'roothtml' --exclude 'CVS' \
		--exclude Make-depend --exclude '*html/gif' \
		--exclude "*tgt_*" --exclude check \
                `cat /tmp/saves` 

htmldocnew:		FORCE
		@for i in $(ALIROOT_DIRS) ; do \
		    echo "Making HTML doc for $$i" ; \
		    rm -rf $(ALICE_ROOT)/$$i/html ; \
		    cd $(ALICE_ROOT)/$$i ; \
		    aliroot -b -q mkhtml.C > /dev/null; \
		    for j in `ls *.C` ; do \
		       echo $$j ; \
		       aliroot -b -q "mkhtml.C(\"$$j\")" > /dev/null; \
		    done \
                done


htmldoc:		FORCE
		@rm -rf html/roothtml
		@rm -f  html/picts
		@rm -f /tmp/macros
		@cd html ;\
		aliroot -q -b "mkhtml.C(0,1)" ;\
		ls ../macros/*.C > /tmp/macros ;\
		for i in $(ALIROOT_DIRS) ; do \
			ls ../$$i/*.C 2>/dev/null >> /tmp/macros ;\
		done ;\
		for i in `cat /tmp/macros` ; do \
			echo $$i ; \
			aliroot -b -q "mkhtml.C(\"$$i\")" > /dev/null ;\
                done ;\
		./makeExampleList ;
		@ln -s ../picts html/picts
		@ln -s ../../picts html/roothtml/picts
		@ln -s ../../../picts html/roothtml/src/picts
		@ln -s ../../../picts html/roothtml/examples/picts

clean:  FORCE
		@rm -f *~ \#*
		@for i in $(ALIROOT_DIRS) ALIROOT ; do \
                    ${MAKE} -C $$i macroclean ; \
                done

libclean:  FORCE
		@rm -f *~ \#*
		@for i in $(LIBRARY_DIRS) ; do \
                    ${MAKE} -C $$i macroclean ; \
                done

allclean: libclean clean

# IRST coding rule check
CHECK_DIRS = $(ALIROOT_DIRS) ALIROOT
check:     
		@for i in $(CHECK_DIRS) ; do \
                    echo "Checking $$i" ; \
                    ${MAKE} -C $$i check ; \
		done

REVENG_DIRS = $(ALIROOT_DIRS)

reveng:
	@for i in $(REVENG_DIRS) ; do \
		echo "Reverse engineering $$i" ; \
		${MAKE} -C $$i reveng ; \
	done








