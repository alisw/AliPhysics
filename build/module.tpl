#-*- Mode: Makefile -*-


ifndef PACKFFLAGS
@PACKAGE@FFLAGS:=$(FFLAGS)
else
@PACKAGE@FFLAGS:=$(PACKFFLAGS)
endif

ifndef PACKCFLAGS
@PACKAGE@CFLAGS:=$(CFLAGS)
else
@PACKAGE@CFLAGS:=$(PACKCFLAGS)
endif

ifndef PACKCXXFLAGS
@PACKAGE@CXXFLAGS:=$(CXXFLAGS)
else
@PACKAGE@CXXFLAGS:=$(PACKCXXFLAGS)
endif

ifndef PACKSOFLAGS
@PACKAGE@SOFLAGS:=$(SOFLAGS)
else
@PACKAGE@SOFLAGS:=$(PACKSOFLAGS)
endif

ifdef DYEXT
ifndef PACKDYFLAGS
@PACKAGE@DYFLAGS:=$(DYFLAGS)
else
@PACKAGE@DYFLAGS:=$(PACKDYFLAGS)
endif
endif

ifndef PACKDCXXFLAGS
@PACKAGE@DCXXFLAGS:=$(CXXFLAGSNO)
else
@PACKAGE@DCXXFLAGS:=$(PACKDCXXFLAGS)
endif


ifdef DHDR
WITHDICT=YES
else
WITHDICT=
endif

# Headerfiles exported to common place:
@PACKAGE@EXPORT:=$(EXPORT)
@PACKAGE@EXPORTDEST:=$(patsubst %,$(EXPORTDIR)/%,$(EXPORT))


#Extra include,libs, defines etc.

@PACKAGE@DEFINE:=$(EDEFINE)

@PACKAGE@INC:=$(patsubst %,-I%,$(EINCLUDE)) -I@MODULE@

@PACKAGE@ELIBS:=$(patsubst %,-l%,$(ELIBS))
@PACKAGE@ELIBSDEP:=$(patsubst %,lib/tgt_$(ALICE_TARGET)/lib%.$(SOEXT),$(ELIBS))
@PACKAGE@ELIBSDIR:=$(patsubst %,-L%,$(ELIBSDIR))

#c sources and headers

@PACKAGE@CS:=$(patsubst %,$(MODDIR)/%,$(CSRCS))
@PACKAGE@CH:=$(patsubst %,$(MODDIR)/%,$(CHDRS))

#Fortran sources 
@PACKAGE@FS:=$(patsubst %,$(MODDIR)/%,$(FSRCS))

#c++ sources and header
@PACKAGE@S:=$(patsubst %,$(MODDIR)/%,$(SRCS))
@PACKAGE@H:=$(patsubst %,$(MODDIR)/%,$(HDRS)) $(EHDRS)

#c++ source subdirectories
@PACKAGE@SDIR:=$(SUBDIR)

#############################################################################
#
#            If special rootcint headerfiles is specified use them
#                         else use all headers

ifndef CINTHDRS
@PACKAGE@CINTHDRS:=$(@PACKAGE@H) 
else
@PACKAGE@CINTHDRS:=$(CINTHDRS)
endif
#############################################################################

# Package Dictionary 

@PACKAGE@DH:=$(MODDIR)/$(DHDR)


#All objects
@PACKAGE@CO:=$(patsubst %,$(MODDIRO)/%, $(CSRCS:.c=.o))
TEMP:=$(FSRCS:.F=.o)
@PACKAGE@FO:=$(patsubst %,$(MODDIRO)/%, $(TEMP:.f=.o))
@PACKAGE@O:= $(patsubst %,$(MODDIRO)/%, $(SRCS:.cxx=.o)) $(@PACKAGE@FO) $(@PACKAGE@CO)



ifdef WITHDICT
  @PACKAGE@DS:=$(MODDIRO)/G__@PACKAGE@.cxx
  @PACKAGE@DO:=$(MODDIRO)/G__@PACKAGE@.o
  @PACKAGE@DDEP:=$(@PACKAGE@DO:.o=.d)
  @PACKAGE@DEP:=$(@PACKAGE@O:.o=.d) $(@PACKAGE@DDEP)
else
  @PACKAGE@DS:=
  @PACKAGE@DO:=
  @PACKAGE@DDEP:=
  @PACKAGE@DEP:=$(@PACKAGE@O:.o=.d)
endif


#The actual library file

@PACKAGE@LIB:=$(LIBPATH)/lib@PACKAGE@.$(SOEXT)

ifneq ($(DYEXT),)
@PACKAGE@DLIB:=$(LIBPATH)/lib@PACKAGE@.$(DYEXT)
endif

@PACKAGE@ALIB:=$(LIBPATH)/lib@PACKAGE@.$(AEXT)

#Add this to the modules libs
@MODULE@LIBS += $(@PACKAGE@LIB)
@MODULE@ALIBS += $(@PACKAGE@ALIB)
ifneq ($(DYEXT),)
@MODULE@DLIBS += $(@PACKAGE@DLIB)
endif

#The actual binary file

@PACKAGE@BIN:=$(BINPATH)/@PACKAGE@

#Add to modules list of binaries
@MODULE@BINS += $(@PACKAGE@BIN)

# Use in the main Makefile

ifeq ($(TYPE),lib)
ALLLIBS += $(@PACKAGE@LIB)
ALLALIBS += $(@PACKAGE@ALIB)
ifneq ($(DYEXT),)
ALLLIBS += $(@PACKAGE@DLIB)
endif
BINLIBS += -l@PACKAGE@
else
ALLEXECS += $(@PACKAGE@BIN)
endif

ifeq ($(DYEXT),)
@PACKAGE@LIB := $(@PACKAGE@LIB)
else
@PACKAGE@LIB := $(@PACKAGE@LIB)
endif

# include all dependence files
INCLUDEFILES +=$(@PACKAGE@DEP)

EXPORTFILES += $(@PACKAGE@EXPORTDEST)

#local rules

#The exportfiles only include if any!!

ifdef @PACKAGE@EXPORT
#$(@PACKAGE@EXPORTDEST): $(patsubst %,@MODULE@/%,$(@PACKAGE@EXPORT))

$(@PACKAGE@EXPORTDEST): $(EXPORTDIR)/%.h: @MODULE@/%.h
ifndef ALIQUIET
	  @echo "***** Copying file $^ to $@ *****"
endif
	  @[ -d $(dir $@) ] || mkdir -p $(dir $@)
	  @cp $^ $@	
endif

#------------------------------------------------------------------------

$(@PACKAGE@LIB):$(@PACKAGE@O) $(@PACKAGE@DO) @MODULE@/module.mk
ifndef ALIQUIET
	  @echo "***** Linking library $@ *****"
endif
	  $(MUTE)TMPDIR=/tmp/@MODULE@$$$$.`date +%M%S` ; \
	  export TMPDIR; mkdir -p $$TMPDIR ; cd $$TMPDIR ; \
	  find $(CURDIR)/@MODULE@/tgt_$(ALICE_TARGET) -name '*.o' -exec ln -s {} . \; ;\
	  \rm -f $(CURDIR)/$@ ;\
	  TMPLIB=$(notdir $(@PACKAGE@LIB)); export TMPLIB;\
	  $(SHLD) $(@PACKAGE@SOFLAGS) -o $(CURDIR)/$@ $(notdir $(@PACKAGE@O) $(@PACKAGE@DO))  $(@PACKAGE@ELIBSDIR) $(@PACKAGE@ELIBS) $(SHLIB);\
	  chmod a-w $(CURDIR)/$@ ;\
	  cd $(ALICE_ROOT) ; \rm -rf $$TMPDIR

ifneq ($(DYEXT),)
$(@PACKAGE@DLIB):$(@PACKAGE@O) $(@PACKAGE@DO) @MODULE@/module.mk
ifndef ALIQUIET
	  @echo "***** Linking library $@ *****"
endif
	  $(MUTE)TMPDIR=/tmp/@MODULE@$$$$.`date +%M%S` ; \
	  export TMPDIR; mkdir -p $$TMPDIR ; cd $$TMPDIR ; \
	  find $(CURDIR)/@MODULE@/tgt_$(ALICE_TARGET) -name '*.o' -exec ln -s {} . \; ;\
	  \rm -f $(CURDIR)/$@ ;\
	  $(DYLD) $(@PACKAGE@DYFLAGS) -o $(CURDIR)/$@ $(notdir $(@PACKAGE@O) $(@PACKAGE@DO))  $(@PACKAGE@ELIBSDIR) $(@PACKAGE@ELIBS) $(DYLIB);\
	  chmod a-w $(CURDIR)/$@ ;\
	  cd $(ALICE_ROOT) ; \rm -rf $$TMPDIR
endif

#------------------------------------------------------------------------

$(@PACKAGE@ALIB):$(@PACKAGE@O) $(@PACKAGE@DO) @MODULE@/module.mk
ifndef ALIQUIET
	  @echo "***** Linking static library $@ *****"
endif
	  $(MUTE)TMPDIR=/tmp/@MODULE@$$$$.`date +%M%S` ; \
	  export TMPDIR; mkdir -p $$TMPDIR ; cd $$TMPDIR ; \
	  find $(CURDIR)/@MODULE@/tgt_$(ALICE_TARGET) -name '*.o' -exec ln -s {} . \; ;\
	  \rm -f $(CURDIR)/$@ ;\
	  TMPLIB=$(notdir $(@PACKAGE@LIB)); export TMPLIB;\
	  $(ALLD) $(ALFLAGS) $(CURDIR)/$@ $(notdir $(@PACKAGE@O) $(@PACKAGE@DO))  $(@PACKAGE@ELIBSDIR) $(@PACKAGE@ELIBS) $(ALLIB);\
      cd $(CURDIR) ; \rm -rf $$TMPDIR
	  $(MUTE)chmod a-w $@


$(@PACKAGE@BIN):$(@PACKAGE@O) $(@PACKAGE@DO) @MODULE@/module.mk $(@PACKAGE@ELIBSDEP)
ifndef ALIQUIET
	  @echo "***** Making executable $@ *****"
endif
ifeq ($(ALIPROFILE),YES)
	$(MUTE)$(LD) $(LDFLAGS) $(@PACKAGE@O) $(ARLIBS) $(SHLIBS) $(LIBS) $(EXEFLAGS) -o $@
else
	  $(MUTE)$(LD) $(LDFLAGS) $(@PACKAGE@O) $(@PACKAGE@DO) $(BINLIBDIRS) $(@PACKAGE@ELIBSDIR) $(@PACKAGE@ELIBS) $(LIBS) $(EXEFLAGS) -o $@
endif

$(@PACKAGE@DS): $(@PACKAGE@CINTHDRS) $(@PACKAGE@DH) @MODULE@/module.mk
ifndef ALIQUIET
	 @echo "***** Creating $@ *****";	
endif
	 @(if [ ! -d '$(dir $@)' ]; then echo "***** Making directory $(dir $@) *****"; mkdir -p $(dir $@); fi;)
	 @\rm -f $(patsubst %.cxx,%.d, $@)
	 $(MUTE)rootcint -f $@ -c $(@PACKAGE@DEFINE) $(CINTFLAGS) $(@PACKAGE@INC) $(@PACKAGE@CINTHDRS) $(@PACKAGE@DH) 

$(@PACKAGE@DO): $(@PACKAGE@DS)
ifndef ALIQUIET
		@echo "***** Compiling $< *****";
endif
		$(MUTE)$(CXX) $(@PACKAGE@DEFINE) -c $(@PACKAGE@INC)  -I$(ALICE_ROOT) $< -o $@ $(@PACKAGE@DCXXFLAGS)

#Different targets for the module

ifeq ($(TYPE),lib)
all-@PACKAGE@: $(@PACKAGE@LIB)
ifneq ($(DYEXT),)
all-@PACKAGE@: $(@PACKAGE@DLIB)
endif
else
all-@PACKAGE@: $(@PACKAGE@BIN)
endif

depend-@PACKAGE@: $(@PACKAGE@DEP)

# determination of object files
$(MODDIRO)/%.o: $(MODDIR)/%.cxx $(MODDIRO)/%.d 
ifndef ALIQUIET
	@echo "***** Compiling $< *****";
endif
	@(if [ ! -d '$(dir $@)' ]; then echo "***** Making directory $(dir $@) *****"; mkdir -p $(dir $@); fi;)
	$(MUTE)$(CXX) $(@PACKAGE@DEFINE) -c $(@PACKAGE@INC)   $< -o $@ $(@PACKAGE@CXXFLAGS)

$(MODDIRO)/%.o: $(MODDIR)/%.F $(MODDIRO)/%.d 
ifndef ALIQUIET
	@echo "***** Compiling $< *****";
endif
	@(if [ ! -d '$(dir $@)' ]; then echo "***** Making directory $(dir $@) *****"; mkdir -p $(dir $@); fi;)
	$(MUTE)$(F77) -c $(@PACKAGE@INC)  $< -o $@ $(@PACKAGE@FFLAGS)

$(MODDIRO)/%.o: $(MODDIR)/%.f $(MODDIRO)/%.d 
ifndef ALIQUIET
	@echo "***** Compiling $< *****";
endif
	@(if [ ! -d '$(dir $@)' ]; then echo "***** Making directory $(dir $@) *****"; mkdir -p $(dir $@); fi;)
	$(MUTE)$(F77) -c $(@PACKAGE@INC)  $< -o $@ $(@PACKAGE@FFLAGS)

$(MODDIRO)/%.o: $(MODDIR)/%.c $(MODDIRO)/%.d 
ifndef ALIQUIET
	@echo "***** Compiling $< *****";
endif
	@(if [ ! -d '$(dir $@)' ]; then echo "***** Making directory $(dir $@) *****"; mkdir -p $(dir $@); fi;)
	$(MUTE)$(CC) $(@PACKAGE@DEFINE) -c  $(@PACKAGE@INC)  $< -o $@   $(@PACKAGE@CFLAGS)

$(@PACKAGE@DDEP): $(@PACKAGE@DS)
ifndef ALIQUIET
		@echo "***** Making dependences for $< *****";
endif
		@(if [ ! -d '$(dir $@)' ]; then echo "***** Making directory $(dir $@) *****"; mkdir -p $(dir $@); fi;)
		@share/alibtool depend "$(@PACKAGE@ELIBSDIR) $(@PACKAGE@INC) $(DEPINC)  $<" > $@

$(MODDIRO)/%.d: $(MODDIRS)/%.cxx
ifndef ALIQUIET
		@echo "***** Making dependences for $< *****";
endif
		@(if [ ! -d '$(dir $@)' ]; then echo "***** Making directory $(dir $@) *****"; mkdir -p $(dir $@); fi;)
		@share/alibtool depend "$(@PACKAGE@DEFINE) $(@PACKAGE@ELIBSDIR) $(@PACKAGE@INC) $(DEPINC)  $<" > $@
$(MODDIRO)/%.d: $(MODDIRS)/%.f
ifndef ALIQUIET
		@echo "***** Making dependences for $< *****";
endif
		@(if [ ! -d '$(dir $@)' ]; then echo "***** Making directory $(dir $@) *****"; mkdir -p $(dir $@); fi;)
		@share/alibtool dependF "$(@PACKAGE@ELIBSDIR) $(@PACKAGE@INC) $(DEPINC)  $<" > $@
$(MODDIRO)/%.d: $(MODDIRS)/%.F
ifndef ALIQUIET
		@echo "***** Making dependences for $< *****";
endif
		@(if [ ! -d '$(dir $@)' ]; then echo "***** Making directory $(dir $@) *****"; mkdir -p $(dir $@); fi;)
		$(MUTE)share/alibtool dependF "$(@PACKAGE@ELIBSDIR) $(@PACKAGE@INC) $(DEPINC)  $<" > $@
$(MODDIRO)/%.d: $(MODDIRS)/%.c
ifndef ALIQUIET
		@echo "***** Making dependences for $< *****";
endif
		@(if [ ! -d '$(dir $@)' ]; then echo "***** Making directory $(dir $@) *****"; mkdir -p $(dir $@); fi;)
		@share/alibtool depend "$(@PACKAGE@DEFINE) $(@PACKAGE@ELIBSDIR) $(@PACKAGE@INC) $(DEPINC) $<" > $@

@PACKAGE@CHECKS := $(patsubst %.cxx,@MODULE@/check/%.viol,$(SRCS))

check-@MODULE@: $(@PACKAGE@CHECKS)

# IRST coding rule check 
@MODULE@/check/%.i : @MODULE@/%.cxx
	@[ -d $(dir $@) ] || mkdir -p $(dir $@)
	$(MUTE)$(CXX) -E $(@PACKAGE@DEFINE) $(@PACKAGE@INC) $< > $@ $(@PACKAGE@CXXFLAGS)
	@cd $(dir $@) ; $(IRST_INSTALLDIR)/patch/patch4alice.prl $(notdir $@)

# IRST coding rule check
@MODULE@/check/$(SUBDIR)/%.viol : @MODULE@/check/$(SUBDIR)/%.i
	@cd @MODULE@ ; [ -r @MODULE@ ] || ln -s ../@MODULE@ @MODULE@
	$(MUTE)echo $@ ; $(CODE_CHECK) $< ./@MODULE@/$(@PACKAGE@SDIR) > $@

# IRST coding rule check
@MODULE@/check/%.viol : @MODULE@/check/%.i
	@cd @MODULE@ ; [ -r @MODULE@ ] || ln -s ../@MODULE@ @MODULE@
	$(MUTE)echo $@ ; $(CODE_CHECK) $< ./@MODULE@ > $@

@PACKAGE@PREPROC       = $(patsubst %.viol,%.i,$(@PACKAGE@CHECKS))

@PACKAGE@REVENGS       = $(patsubst %.viol,%.ii,$(@PACKAGE@CHECKS))

.SECONDARY: $(@PACKAGE@REVENGS) $(@PACKAGE@PREPROC)

reveng-@PACKAGE@:		@PACKAGE@/check/classDiagram.dot

@PACKAGE@/check/classDiagram.dot:	$(@PACKAGE@PREPROC)
	@$(REV_ENG) $^
	@-mv classDiagram.dot $@

revdisp-@PACKAGE@:	reveng-@PACKAGE@
	@echo revdisp for @PACKAGE@
	@cd @PACKAGE@/check ; \
      $(IRST_INSTALLDIR)/webreveng/create-class-diagram-pages.sh
	@sed -e "s/\@PACKAGE\@/@PACKAGE@/g" < $(ALICE_ROOT)/build/HomePage.html > @PACKAGE@/check/HomePage.html

