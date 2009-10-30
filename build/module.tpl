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

ifndef PACKLDFLAGS
@PACKAGE@LDFLAGS:=$(LDFLAGS)
else
@PACKAGE@LDFLAGS:=$(PACKLDFLAGS)
endif

ifndef PACKDCXXFLAGS
@PACKAGE@DCXXFLAGS:=$(CXXFLAGSNO)
else
@PACKAGE@DCXXFLAGS:=$(PACKDCXXFLAGS)
endif

ifndef PACKBLIBS
@PACKAGE@BLIBS:=$(LIBS)
else
@PACKAGE@BLIBS:=$(PACKBLIBS)
endif

WITHDICT=
ifneq ($(DHDR),)
WITHDICT=YES
else
ifneq ($(CINTAUTOLINK),)
WITHDICT=YES
endif
endif

# Headerfiles exported to common place:
@PACKAGE@EXPORT:=$(EXPORT)
@PACKAGE@EXPORTDEST:=$(patsubst %,$(EXPORTDIR)/%,$(EXPORT))


#Extra include,libs, defines etc.

@PACKAGE@DEFINE:=$(EDEFINE) -D_MODULE_=\"@MODULE@\"

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

#############################################################################
#
#            If special rootcint headerfiles is specified use them
#                         else use all headers

ifndef CINTHDRS
@PACKAGE@CINTHDRS:=$(@PACKAGE@H)
else
@PACKAGE@CINTHDRS:=$(CINTHDRS:%=@MODULE@/%)
#@PACKAGE@CINTHDRS:=$(pathsubst %,@MODULE@/%,$(CINTHDRS))
endif
@PACKAGE@CINTCLASSES:=$(patsubst %.h,%,$(notdir $(@PACKAGE@CINTHDRS)))
#############################################################################

# Package Dictionary 
ifneq ($(DHDR),)
@PACKAGE@DH:=$(MODDIR)/$(DHDR)
else
@PACKAGE@DH:=
endif
ifneq ($(CINTAUTOLINK),)
@PACKAGE@DAL:=$(MODDIRO)/G__@PACKAGE@AutoLinkDef.h
@PACKAGE@DH+=$(@PACKAGE@DAL)
endif

#All objects
@PACKAGE@CXXO:=$(patsubst %,$(MODDIRO)/%, $(SRCS:.cxx=.o))
@PACKAGE@CO:=$(patsubst %,$(MODDIRO)/%, $(CSRCS:.c=.o))
@PACKAGE@SMALLFO:=$(patsubst %.f,$(MODDIRO)/%.o, $(filter %.f, $(FSRCS)))
@PACKAGE@CAPITFO:=$(patsubst %.F,$(MODDIRO)/%.o, $(filter %.F, $(FSRCS)))
@PACKAGE@FO:=$(@PACKAGE@SMALLFO) $(@PACKAGE@CAPITFO)
@PACKAGE@O:= $(@PACKAGE@CXXO) $(@PACKAGE@FO) $(@PACKAGE@CO)
@PACKAGE@CHV:=$(patsubst %,$(MODDIRC)/%, $(SRCS:.cxx=.viol))
@PACKAGE@SML:=$(sort $(patsubst %,$(MODDIRZ)/%, $(filter $(SRCS:.cxx=.smell),$(HDRS:.h=.smell))))


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

@PACKAGE@ALIB:=$(LIBPATH)/lib@PACKAGE@.$(AEXT)

#Add this to the modules libs
ifeq ($(TYPE),lib)
@MODULE@LIBS += $(@PACKAGE@LIB)
@MODULE@ALIBS += $(@PACKAGE@ALIB)
endif

#The actual binary file

@PACKAGE@BIN:=$(BINPATH)/@PACKAGE@

#Add to modules list of binaries
@MODULE@BINS += $(@PACKAGE@BIN)

# Use in the main Makefile

ifeq ($(TYPE),lib)
ALLLIBS += $(@PACKAGE@LIB)
ALLALIBS += $(@PACKAGE@ALIB)
BINLIBS += -l@PACKAGE@
else
ALLEXECS += $(@PACKAGE@BIN)
endif

@PACKAGE@LIB := $(@PACKAGE@LIB)

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
	  \rm -f "$(CURDIR)"/$@ ;\
	  cd $(@MODULE@DIRO) ;\
	  $(SHLD) $(@PACKAGE@SOFLAGS) -o "$(CURDIR)"/$@ $(patsubst $(@MODULE@DIRO)/%,%,$(@PACKAGE@O) $(@PACKAGE@DO))  $(@PACKAGE@ELIBSDIR) $(@PACKAGE@ELIBS) $(SHLIB);\
	  chmod a-w "$(CURDIR)"/$@ ;\
	  cd $(ALICE_ROOT)

#------------------------------------------------------------------------

$(@PACKAGE@ALIB):$(@PACKAGE@O) $(@PACKAGE@DO) @MODULE@/module.mk
ifndef ALIQUIET
	  @echo "***** Linking static library $@ *****"
endif
	  \rm -f "$(CURDIR)"/$@ ;\
	  cd $(@MODULE@DIRO) ;\
	  $(ALLD) $(ALFLAGS) "$(CURDIR)"/$@ $(patsubst $(@MODULE@DIRO)/%,%,$(@PACKAGE@O) $(@PACKAGE@DO))  $(@PACKAGE@ELIBSDIR) $(@PACKAGE@ELIBS) $(ALLIB);\
          cd $(ALICE_ROOT)


$(@PACKAGE@BIN):$(@PACKAGE@O) $(@PACKAGE@DO) @MODULE@/module.mk
ifndef ALIQUIET
	  @echo "***** Making executable $@ *****"
endif
ifeq ($(ALIPROFILE),YES)
	$(MUTE)$(LD) $(@PACKAGE@LDFLAGS) $(@PACKAGE@O) $(ARLIBS) $(SHLIBS) $(@PACKAGE@BLIBS) $(EXEFLAGS) -o $@
else
	  $(MUTE)$(LD) $(@PACKAGE@LDFLAGS) $(@PACKAGE@O) $(@PACKAGE@DO) $(BINLIBDIRS) $(@PACKAGE@ELIBSDIR) $(@PACKAGE@ELIBS) $(@PACKAGE@BLIBS) $(EXEFLAGS) -o $@
endif

$(@PACKAGE@DAL): $(@PACKAGE@CINTHDRS) @MODULE@/module.mk
ifndef ALIQUIET
	 @echo "***** Creating $@ *****";
endif
	$(MUTE)echo '//automatically generated ROOT DICT definition' > $@
	$(MUTE)echo '//!!! DO NOT EDIT THIS FILE !!!' >> $@
	$(MUTE)echo '#ifdef __CINT__' >> $@
	$(MUTE)echo '#pragma link off all globals;' >> $@
	$(MUTE)echo '#pragma link off all classes;' >> $@
	$(MUTE)echo '#pragma link off all functions;' >> $@
	$(MUTE)$(foreach i, $(@PACKAGE@CINTCLASSES), \
	   echo "#pragma link C++ class $(i)+;" >> $@ ;)
	$(MUTE)echo '#endif' >> $@

ifneq ($(findstring solaris,$(ALICE_TARGET)),)
# shameful hack....
CINTFLAGS := $(filter-out -library=stlport4,$(CINTFLAGS))
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
all-@MODULE@: $(@PACKAGE@LIB)
else
all-@MODULE@: $(@PACKAGE@BIN)
endif


depend-@PACKAGE@: $(@PACKAGE@DEP)

# determination of object files
$(@PACKAGE@CXXO): $(MODDIRO)/%.o: $(MODDIR)/%.cxx $(MODDIRO)/%.d 
ifndef ALIQUIET
	@echo "***** Compiling $< *****";
endif
	@(if [ ! -d '$(dir $@)' ]; then echo "***** Making directory $(dir $@) *****"; mkdir -p $(dir $@); fi;)
	$(MUTE)$(CXX) $(@PACKAGE@DEFINE) -c $(@PACKAGE@INC)   $< -o $@ $(@PACKAGE@CXXFLAGS)

$(@PACKAGE@CAPITFO): $(MODDIRO)/%.o: $(MODDIR)/%.F $(MODDIRO)/%.d 
ifndef ALIQUIET
	@echo "***** Compiling $< *****";
endif
	@(if [ ! -d '$(dir $@)' ]; then echo "***** Making directory $(dir $@) *****"; mkdir -p $(dir $@); fi;)
	$(MUTE)$(F77) -c $(@PACKAGE@INC)  $< -o $@ $(@PACKAGE@FFLAGS)

$(@PACKAGE@SMALLFO): $(MODDIRO)/%.o: $(MODDIR)/%.f $(MODDIRO)/%.d 
ifndef ALIQUIET
	@echo "***** Compiling $< *****";
endif
	@(if [ ! -d '$(dir $@)' ]; then echo "***** Making directory $(dir $@) *****"; mkdir -p $(dir $@); fi;)
	$(MUTE)$(F77) -c $(@PACKAGE@INC)  $< -o $@ $(@PACKAGE@FFLAGS)

$(@PACKAGE@CO): $(MODDIRO)/%.o: $(MODDIR)/%.c $(MODDIRO)/%.d 
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

$(@PACKAGE@CXXO:.o=.d): $(MODDIRO)/%.d: $(MODDIRS)/%.cxx
ifndef ALIQUIET
		@echo "***** Making dependences for $< *****";
endif
		@(if [ ! -d '$(dir $@)' ]; then echo "***** Making directory $(dir $@) *****"; mkdir -p $(dir $@); fi;)
		@share/alibtool depend "$(@PACKAGE@DEFINE) $(@PACKAGE@ELIBSDIR) $(@PACKAGE@INC) $(DEPINC)  $<" > $@

$(@PACKAGE@SMALLFO:.o=.d): $(MODDIRO)/%.d: $(MODDIRS)/%.f
ifndef ALIQUIET
		@echo "***** Making dependences for $< *****";
endif
		@(if [ ! -d '$(dir $@)' ]; then echo "***** Making directory $(dir $@) *****"; mkdir -p $(dir $@); fi;)
		@share/alibtool dependF "$(@PACKAGE@ELIBSDIR) $(@PACKAGE@INC) $(DEPINC)  $<" > $@

$(@PACKAGE@CAPITFO:.o=.d): $(MODDIRO)/%.d: $(MODDIRS)/%.F
ifndef ALIQUIET
		@echo "***** Making dependences for $< *****";
endif
		@(if [ ! -d '$(dir $@)' ]; then echo "***** Making directory $(dir $@) *****"; mkdir -p $(dir $@); fi;)
		$(MUTE)share/alibtool dependF "$(@PACKAGE@ELIBSDIR) $(@PACKAGE@INC) $(DEPINC)  $<" > $@

$(@PACKAGE@CO:.o=.d): $(MODDIRO)/%.d: $(MODDIRS)/%.c
ifndef ALIQUIET
		@echo "***** Making dependences for $< *****";
endif
		@(if [ ! -d '$(dir $@)' ]; then echo "***** Making directory $(dir $@) *****"; mkdir -p $(dir $@); fi;)
		@share/alibtool depend "$(@PACKAGE@DEFINE) $(@PACKAGE@ELIBSDIR) $(@PACKAGE@INC) $(DEPINC) $<" > $@

.PRECIOUS: $(patsubst %.cxx,$(MODDIRO)/%.d,$(SRCS))
.PRECIOUS: $(patsubst %.c,$(MODDIRO)/%.d,$(CSRCS))
.PRECIOUS: $(patsubst %.F,$(MODDIRO)/%.d,$(patsubst %.f,$(MODDIRO)/%.d,$(FSRCS)))

check-@MODULE@: $(@PACKAGE@CHV)

# IRST coding rule check 
$(@PACKAGE@CHV:.viol=.i): $(MODDIRC)/%.i: $(MODDIR)/%.cxx $(MODDIRO)/%.d
	@[ -d $(dir $@) ] || mkdir -p $(dir $@)
	$(MUTE)$(CXX) -E $(@PACKAGE@DEFINE) $(@PACKAGE@INC) -I. $< > $@ $(@PACKAGE@CXXFLAGS)
	@cd $(dir $@) ; $(IRST_INSTALLDIR)/patch/patch4alice.prl $(notdir $@)

# IRST coding rule check
$(@PACKAGE@CHV): $(MODDIRC)/%.viol: $(MODDIRC)/%.i
	$(MUTE)echo $@ ; $(CODE_CHECK) $< $(shell echo $(dir $<) | sed -e 's:/check::') > $@

.SECONDARY: $(@PACKAGE@CHV:.viol=.i) $(@PACKAGE@CHI:.viol=.ii)

PACKREVENG += $(@PACKAGE@CHV:.viol=.ii)

# IRST code smell checker

smell-@MODULE@: $(@PACKAGE@SML)

$(@PACKAGE@SML:.smell=_cxx.ml) : $(MODDIRZ)/%_cxx.ml : $(MODDIR)/%.cxx
	@[ -d $(dir $@) ] || mkdir -p $(dir $@)
	$(MUTE)src2srcml $< $@

$(@PACKAGE@SML:.smell=_h.ml) : $(MODDIRZ)/%_h.ml : $(MODDIR)/%.h
	@[ -d $(dir $@) ] || mkdir -p $(dir $@)
	$(MUTE)src2srcml $< $@

$(@PACKAGE@SML) : $(MODDIRZ)/%.smell : $(MODDIRZ)/%_cxx.ml $(MODDIRZ)/%_h.ml
	$(MUTE)echo smelling $@
	$(MUTE)java -classpath $(SMELL_DETECTOR_DIR):$(SMELL_DETECTOR_DIR)/xom-1.1.jar -Xmx500m SmellDetector $^ > $@
	$(MUTE)[ -s $@ ] || touch $@

.SECONDARY: $(@PACKAGE@SML:.smell=_cxx.ml) $(@PACKAGE@SML:.smell=_h.ml)

# targets to create .par archives (jgrosseo)
@PACKAGE@.par: $(patsubst %,@MODULE@/@PACKAGE@/%,$(filter-out dict.%, $(HDRS) $(SRCS) $(DHDR) $(PKGFILE) $(FSRCS) Makefile Makefile.arch lib@PACKAGE@.pkg PROOF-INF))
	@echo "Creating archive" $@ ...
	@cd @MODULE@; (tar cfzh ../$@ @PACKAGE@ 2> /dev/null && echo "package" $@ "created in" $(PWD)/$@) || (tar cfzh /tmp/$@ @PACKAGE@ 2> /dev/null && echo "package" $@ "created in /tmp/"$@)
	@rm -rf @MODULE@/@PACKAGE@

@MODULE@/@PACKAGE@/Makefile: @MODULE@/Makefile
	@echo Copying $< to $@ with transformations
	@[ -d $(dir $@) ] || mkdir -p $(dir $@)
	@sed 's/include \$$(ROOTSYS)\/test\/Makefile.arch/include Makefile.arch/; s/PACKAGE = .*/PACKAGE = @PACKAGE@/' < $^ > $@

@MODULE@/@PACKAGE@/Makefile.arch: $(ROOTSYS)/test/Makefile.arch
	@echo Copying $< to $@
	@[ -d $(dir $@) ] || mkdir -p $(dir $@)
	@cp -pR $^ $@

@MODULE@/@PACKAGE@/PROOF-INF: @MODULE@/PROOF-INF.@PACKAGE@
	@echo Copying $< to $@
	@[ -d $(dir $@) ] || mkdir -p $(dir $@)
	@cp -pR $^ $@

@MODULE@/@PACKAGE@/%: @MODULE@/%
	@echo Copying $< to $@
	@[ -d $(dir $@) ] || mkdir -p $(dir $@)
	@cp -pR $< $@

test-@PACKAGE@.par: @PACKAGE@.par
	@echo "INFO: The file $< is now tested, in case of an error check in par-tmp/@PACKAGE@."
	@mkdir -p par-tmp
	@cd par-tmp; tar xfz ../$<;	cd $(subst .par,,$<); PROOF-INF/BUILD.sh
	@rm -rf par-tmp/@PACKAGE@
	@echo "INFO: Testing succeeded (already cleaned up)"
