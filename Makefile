############################### Main Makefile #################################

# Include machine specific definitions

include $(ALICE_ROOT)/conf/GeneralDef
include $(ALICE_ROOT)/conf/MachineDef.$(ALICE_TARGET)

MAKEFLAGS =

##### MACROS #####

PACKAGE = Main

##### Module libraries #####

ALIROOT_DIRS		= STEER TGeant3 TRD PHOS TPC ZDC MUON PMD FMD TOF ITS CASTOR \
		  RICH STRUCT EVGEN

##### TARGETS #####
 
default:      lib bin alilibs aliroot

lib bin:
	@mkdir $@

alilibs:  lib
	@for i in $(ALIROOT_DIRS) ; do \
	   ${MAKE} -C $$i depend ; \
	done
	@for i in $(ALIROOT_DIRS) ; do \
	   ${MAKE} -C $$i ; \
	done

aliroot: bin
	@${MAKE} -C ALIROOT

geant321:  lib
	@-${MAKE} -C GEANT321 depend
	@${MAKE} -C GEANT321

pythia:    lib
	@-${MAKE} -C PYTHIA depend
	@${MAKE} -C PYTHIA

pdf:       lib
	@-${MAKE} -C PDF depend
	@${MAKE} -C PDF

minicern:  lib
	@-${MAKE} -C MINICERN depend
	@${MAKE} -C MINICERN

cernlibs: geant321 pythia minicern pdf

all:	cernlibs default

FORCE:

############################### General Macros ################################

include $(ALICE_ROOT)/conf/GeneralMacros

############################### Specific Macros ###############################

STRUCT_DIRS	= html conf macros data share include Euclid picts \
                  Makefile README .rootrc

LIBRARY_DIRS	= MINICERN GEANT321 PYTHIA PDF

dist: AliRoot3.01.tar.gz

AliRoot3.01.tar.gz: $(STRUCT_DIRS) $(ALIROOT_DIRS) ALIROOT

distall: AliOffline3.01.tar.gz

AliOffline3.01.tar.gz: $(STRUCT_DIRS) $(ALIROOT_DIRS) $(LIBRARY_DIRS) ALIROOT

distlib: AliLibs3.01.tar.gz

AliLibs3.01.tar.gz: $(LIBRARY_DIRS)

AliRoot3.01.tar.gz AliLibs3.01.tar.gz AliOffline3.01.tar.gz:
		@rm -f $(ALICE)/$@ 
		@rm -f `find . -name '*~' -print` \
                       `find . -name '*.bak' -print` \
                       `find . -name '.*~' -print` \
		       `find . -name '*\#*' -print` 
		@rm -f /tmp/saves
		@ls -1d $^ | sed -e "s/^/$(ALICE_LEVEL)\//" > /tmp/saves
	 	@cd $(ALICE) ; \
                gtar cvfz $@ --exclude '*.o' --exclude '*Cint.*' \
                --exclude '*/roothtml/*' --exclude '*/CVS' \
		--exclude Make-depend --exclude '*html/gif' \
                `cat /tmp/saves` 

alidepend:
		@for i in $(ALIROOT_DIRS) ; do \
                    ${MAKE} -C $$i depend ; \
                done

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
