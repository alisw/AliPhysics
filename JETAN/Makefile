# $Id$
#
# Makefile for the jetan package (standalone)

ifdef ALIVERBOSE
MUTE:=
else
MUTE:=@
endif

#-------------------------------------------------------------------------------
# Include library definition

include $(ALICE_ROOT)/build/Makefile.config

#-------------------------------------------------------------------------------
# Location where to install libraries and binaries and common header files

LIBPATH      := ./lib/tgt_$(ALICE_TARGET)
BINPATH      := ./bin/tgt_$(ALICE_TARGET)
EXPORTDIR    := $(ALICE_ROOT)/include
BINLIBDIRS   := -L$(ALICE_ROOT)/$(LIBPATH)

#-------------------------------------------------------------------------------
# Include machine dependent macros

include $(ALICE_ROOT)/build/Makefile.$(ALICE_TARGET)

#-------------------------------------------------------------------------------
# Check if called with debug

ifeq ($(ALIDEBUG),YES)
override ALICE_TARGET := $(ALICE_TARGET)DEBUG
FFLAGS   := -g $(filter-out -O%,$(FFLAGS))
CXXFLAGS := -g $(filter-out -O%,$(CXXLAGS))
CFLAGS   := -g $(filter-out -O%,$(CLAGS))
SOFLAGS  := -g $(filter-out -O%,$(SOFLAGS))
LDFLAGS  := -g $(filter-out -O%,$(LDFLAGS))
endif

#-------------------------------------------------------------------------------
# Check if called with profile

ifeq ($(ALIPROFILE),YES)
override ALICE_TARGET := $(ALICE_TARGET)PROF
FFLAGS   += -pg
CXXFLAGS += -pg
CFLAGS   += -pg
SOFLAGS  += -pg
LDFLAGS  += -pg
endif

#-------------------------------------------------------------------------------
# Check if DATE is installed

ifdef DATE_ROOT
DATEFLAGS  = -DALI_DATE -D${DATE_SYS} -DDATE_SYS=${DATE_SYS} -Dlong32=${DATE_LONG32} \
             -Dlong64=${DATE_LONG64} -DdatePointer=${DATE_POINTER}
CXXFLAGS  += $(DATEFLAGS)
CFLAGS    += $(DATEFLAGS)
CINTFLAGS += $(DATEFLAGS)
DEPINC    += $(DATEFLAGS)
endif

#-------------------------------------------------------------------------------
# Default include dirs for C++, Fortran, Cint, and dependencies
# The module directory will be added by each module

GENINC     := -I$(ALICE_ROOT)/include -I$(shell root-config --incdir)
CXXFLAGS   += $(GENINC)
CXXFLAGSNO += $(GENINC)
CINTFLAGS  += $(GENINC)
DEPINC     += $(GENINC)

LIBS := $(GLIBS) $(ROOTLIBS) $(SYSLIBS)

#Determine if it's a library or a executable
TYPE=lib

# Package head directory, source and include directories
CURDIR=$(shell pwd)
MODDIR:=.
MODDIRS:=$(MODDIR)
MODDIRI:=$(MODDIR)
MODDIRO:=./.objs

# Reseting variables before importing pkg-file
SRCS:=
HDRS:=
FSRCS:=
DHDR:=
CSRCS:=
CHDRS:=
EINCLUDE:=
EDEFINE:=
ELIBS:=
ELIBSDIR:=
PACKFFLAGS:=
PACKCXXFLAGS:=
PACKCFLAGS:=
PACKDYFLAGS:=
PACKSOFLAGS:=
EXPORT:=
EHDRS:=
CINTHDRS:=
ARLIBS:=
SHLIBS:=
SUBDIR:=

include libJETAN.pkg

ifndef PACKFFLAGS
JETANFFLAGS:=$(FFLAGS)
else
JETANFFLAGS:=$(PACKFFLAGS)
endif

ifndef PACKCFLAGS
JETANCFLAGS:=$(CFLAGS)
else
JETANCFLAGS:=$(PACKCFLAGS)
endif

ifndef PACKCXXFLAGS
JETANCXXFLAGS:=$(CXXFLAGS)
else
JETANCXXFLAGS:=$(PACKCXXFLAGS)
endif

ifndef PACKSOFLAGS
JETANSOFLAGS:=$(SOFLAGS)
else
JETANSOFLAGS:=$(PACKSOFLAGS)
endif

ifdef DYEXT
ifndef PACKDYFLAGS
JETANDYFLAGS:=$(DYFLAGS)
else
JETANDYFLAGS:=$(PACKDYFLAGS)
endif
endif

ifndef PACKDCXXFLAGS
JETANDCXXFLAGS:=$(CXXFLAGSNO)
else
JETANDCXXFLAGS:=$(PACKDCXXFLAGS)
endif


ifdef DHDR
WITHDICT=YES
else
WITHDICT=
endif

# Headerfiles exported to common place:
JETANEXPORT:=$(EXPORT)
JETANEXPORTDEST:=$(patsubst %,$(EXPORTDIR)/%,$(EXPORT))


#Extra include,libs, defines etc.

JETANDEFINE:=$(EDEFINE)

JETANINC:=$(patsubst %,-I%,$(EINCLUDE)) -I.

JETANELIBS:=$(patsubst %,-l%,$(ELIBS))
JETANELIBSDEP:=$(patsubst %,lib/tgt_$(ALICE_TARGET)/lib%.$(SOEXT),$(ELIBS))
JETANELIBSDIR:=$(patsubst %,-L%,$(ELIBSDIR))

#c sources and headers

JETANCS:=$(patsubst %,$(MODDIR)/%,$(CSRCS))
JETANCH:=$(patsubst %,$(MODDIR)/%,$(CHDRS))

#Fortran sources 
JETANFS:=$(patsubst %,$(MODDIR)/%,$(FSRCS))

#c++ sources and header
JETANS:=$(patsubst %,$(MODDIR)/%,$(SRCS))
JETANH:=$(patsubst %,$(MODDIR)/%,$(HDRS)) $(EHDRS)

#c++ source subdirectories
JETANSDIR:=$(SUBDIR)

#############################################################################
#
#            If special rootcint headerfiles is specified use them
#                         else use all headers

ifndef CINTHDRS
JETANCINTHDRS:=$(JETANH) 
else
JETANCINTHDRS:=$(CINTHDRS)
endif
#############################################################################

# Package Dictionary 

JETANDH:=$(MODDIR)/$(DHDR)


#All objects
JETANCO:=$(patsubst %,$(MODDIRO)/%, $(CSRCS:.c=.o))
TEMP:=$(FSRCS:.F=.o)
JETANFO:=$(patsubst %,$(MODDIRO)/%, $(TEMP:.f=.o))
JETANO:= $(patsubst %,$(MODDIRO)/%, $(SRCS:.cxx=.o)) $(JETANFO) $(JETANCO)



ifdef WITHDICT
  JETANDS:=$(MODDIRO)/G__JETAN.cxx
  JETANDO:=$(MODDIRO)/G__JETAN.o
  JETANDDEP:=$(JETANDO:.o=.d)
  JETANDEP:=$(JETANO:.o=.d) $(JETANDDEP)
else
  JETANDS:=
  JETANDO:=
  JETANDDEP:=
  JETANDEP:=$(JETANO:.o=.d)
endif


#The actual library file

JETANLIB:=$(LIBPATH)/libJETAN.$(SOEXT)

ifneq ($(DYEXT),)
JETANDLIB:=$(LIBPATH)/libJETAN.$(DYEXT)
endif

JETANALIB:=$(LIBPATH)/libJETAN.$(AEXT)

#Add this to the modules libs
JETANLIBS += $(JETANLIB)
JETANALIBS += $(JETANALIB)
ifneq ($(DYEXT),)
JETANDLIBS += $(JETANDLIB)
endif

#The actual binary file

JETANBIN:=$(BINPATH)/JETAN

#Add to modules list of binaries
JETANBINS += $(JETANBIN)

# Use in the main Makefile

ifeq ($(TYPE),lib)
ALLLIBS += $(JETANLIB)
ALLALIBS += $(JETANALIB)
ifneq ($(DYEXT),)
ALLLIBS += $(JETANDLIB)
endif
BINLIBS += -lJETAN
else
ALLEXECS += $(JETANBIN)
endif

ifeq ($(DYEXT),)
JETANLIB := $(JETANLIB)
else
JETANLIB := $(JETANLIB)
endif

# include all dependence files
INCLUDEFILES +=$(JETANDEP)

EXPORTFILES += $(JETANEXPORTDEST)

#local rules

#The exportfiles only include if any!!

ifdef JETANEXPORT
#$(JETANEXPORTDEST): $(patsubst %,JETAN/%,$(JETANEXPORT))

$(JETANEXPORTDEST): $(EXPORTDIR)/%.h: JETAN/%.h
ifndef ALIQUIET
	  @echo "***** Copying file $^ to $@ *****"
endif
	  @[ -d $(dir $@) ] || mkdir -p $(dir $@)
	  @cp $^ $@	
endif

#------------------------------------------------------------------------

$(JETANLIB):$(JETANO) $(JETANDO)
ifndef ALIQUIET
	  @echo "***** Linking library $@ *****"
endif
	$(MUTE)mkdir -p $(LIBPATH)
	$(MUTE)TMPDIR=/tmp/JETAN$$$$.`date +%M%S` ; \
	export TMPDIR; mkdir -p $$TMPDIR ; cd $$TMPDIR ; \
	find $(CURDIR)/.objs/ -name '*.o' -exec ln -s {} . \; ;\
	\rm -f $(CURDIR)/$@ ;\
	TMPLIB=$(notdir $(JETANLIB)); export TMPLIB;\
	$(SHLD) $(JETANSOFLAGS) -o $(CURDIR)/$@ $(notdir $(JETANO) $(JETANDO))  $(JETANELIBSDIR) $(JETANELIBS) $(SHLIB);\
	chmod a-w $(CURDIR)/$@ ;\
	cd $(ALICE_ROOT) ; \rm -rf $$TMPDIR

ifneq ($(DYEXT),)
$(JETANDLIB):$(JETANO) $(JETANDO)
ifndef ALIQUIET
	  @echo "***** Linking library $@ *****"
endif
	$(MUTE)mkdir -p $(LIBPATH)
	$(MUTE)TMPDIR=/tmp/JETAN$$$$.`date +%M%S` ; \
	export TMPDIR; mkdir -p $$TMPDIR ; cd $$TMPDIR ; \
	find $(CURDIR)/.objs/ -name '*.o' -exec ln -s {} . \; ;\
	\rm -f $(CURDIR)/$@ ;\
	$(DYLD) $(JETANDYFLAGS) -o $(CURDIR)/$@ $(notdir $(JETANO) $(JETANDO))  $(JETANELIBSDIR) $(JETANELIBS) $(DYLIB);\
	chmod a-w $(CURDIR)/$@ ;\
	cd $(ALICE_ROOT) ; \rm -rf $$TMPDIR
endif

#------------------------------------------------------------------------

$(JETANALIB):$(JETANO) $(JETANDO)
ifndef ALIQUIET
	@echo "***** Linking static library $@ *****"
endif
	$(MUTE)mkdir -p $(LIBPATH)
	$(MUTE)TMPDIR=/tmp/JETAN$$$$.`date +%M%S` ; \
	export TMPDIR; mkdir -p $$TMPDIR ; cd $$TMPDIR ; \
	find $(CURDIR)/JETAN/tgt_$(ALICE_TARGET) -name '*.o' -exec ln -s {} . \; ;\
	\rm -f $(CURDIR)/$@ ;\
	TMPLIB=$(notdir $(JETANLIB)); export TMPLIB;\
	$(ALLD) $(ALFLAGS) $(CURDIR)/$@ $(notdir $(JETANO) $(JETANDO))  $(JETANELIBSDIR) $(JETANELIBS) $(ALLIB);\
	cd $(CURDIR) ; \rm -rf $$TMPDIR
	$(MUTE)chmod a-w $@


$(JETANBIN):$(JETANO) $(JETANDO) $(JETANELIBSDEP)
ifndef ALIQUIET
	@echo "***** Making executable $@ *****"
endif
ifeq ($(ALIPROFILE),YES)
	$(MUTE)$(LD) $(LDFLAGS) $(JETANO) $(ARLIBS) $(SHLIBS) $(LIBS) $(EXEFLAGS) -o $@
else
	$(MUTE)$(LD) $(LDFLAGS) $(JETANO) $(JETANDO) $(BINLIBDIRS) $(JETANELIBSDIR) $(JETANELIBS) $(LIBS) $(EXEFLAGS) -o $@
endif

$(JETANDS): $(JETANCINTHDRS) $(JETANDH)
ifndef ALIQUIET
	@echo "***** Creating $@ *****";	
endif
	@(if [ ! -d '$(dir $@)' ]; then echo "***** Making directory $(dir $@) *****"; mkdir -p $(dir $@); fi;)
	@\rm -f $(patsubst %.cxx,%.d, $@)
	$(MUTE)rootcint -f $@ -c $(JETANDEFINE) $(CINTFLAGS) $(JETANINC) $(JETANCINTHDRS) $(JETANDH) 

$(JETANDO): $(JETANDS)
ifndef ALIQUIET
	@echo "***** Compiling $< *****";
endif
	$(MUTE)$(CXX) $(JETANDEFINE) -c $(JETANINC)  -I$(ALICE_ROOT) $< -o $@ $(JETANDCXXFLAGS)

#Different targets for the module

ifeq ($(TYPE),lib)
all-JETAN: $(JETANLIB)
ifneq ($(DYEXT),)
all-JETAN: $(JETANDLIB)
endif
else
all-JETAN: $(JETANBIN)
endif

depend-JETAN: $(JETANDEP)

# determination of object files
$(MODDIRO)/%.o: $(MODDIR)/%.cxx $(MODDIRO)/%.d 
ifndef ALIQUIET
	@echo "***** Compiling $< *****";
endif
	@(if [ ! -d '$(dir $@)' ]; then echo "***** Making directory $(dir $@) *****"; mkdir -p $(dir $@); fi;)
	$(MUTE)$(CXX) $(JETANDEFINE) -c $(JETANINC)   $< -o $@ $(JETANCXXFLAGS)

$(MODDIRO)/%.o: $(MODDIR)/%.F $(MODDIRO)/%.d 
ifndef ALIQUIET
	@echo "***** Compiling $< *****";
endif
	@(if [ ! -d '$(dir $@)' ]; then echo "***** Making directory $(dir $@) *****"; mkdir -p $(dir $@); fi;)
	$(MUTE)$(F77) -c $(JETANINC)  $< -o $@ $(JETANFFLAGS)

$(MODDIRO)/%.o: $(MODDIR)/%.f $(MODDIRO)/%.d 
ifndef ALIQUIET
	@echo "***** Compiling $< *****";
endif
	@(if [ ! -d '$(dir $@)' ]; then echo "***** Making directory $(dir $@) *****"; mkdir -p $(dir $@); fi;)
	$(MUTE)$(F77) -c $(JETANINC)  $< -o $@ $(JETANFFLAGS)

$(MODDIRO)/%.o: $(MODDIR)/%.c $(MODDIRO)/%.d 
ifndef ALIQUIET
	@echo "***** Compiling $< *****";
endif
	@(if [ ! -d '$(dir $@)' ]; then echo "***** Making directory $(dir $@) *****"; mkdir -p $(dir $@); fi;)
	$(MUTE)$(CC) $(JETANDEFINE) -c  $(JETANINC)  $< -o $@   $(JETANCFLAGS)

$(JETANDDEP): $(JETANDS)
ifndef ALIQUIET
		@echo "***** Making dependences for $< *****";
endif
		@(if [ ! -d '$(dir $@)' ]; then echo "***** Making directory $(dir $@) *****"; mkdir -p $(dir $@); fi;)
		@$(ALICE_ROOT)/share/alibtool depend "$(JETANELIBSDIR) $(JETANINC) $(DEPINC)  $<" > $@

$(MODDIRO)/%.d: $(MODDIRS)/%.cxx
ifndef ALIQUIET
		@echo "***** Making dependences for $< *****";
endif
		@(if [ ! -d '$(dir $@)' ]; then echo "***** Making directory $(dir $@) *****"; mkdir -p $(dir $@); fi;)
		@$(ALICE_ROOT)/share/alibtool depend "$(JETANDEFINE) $(JETANELIBSDIR) $(JETANINC) $(DEPINC)  $<" > $@
$(MODDIRO)/%.d: $(MODDIRS)/%.f
ifndef ALIQUIET
		@echo "***** Making dependences for $< *****";
endif
		@(if [ ! -d '$(dir $@)' ]; then echo "***** Making directory $(dir $@) *****"; mkdir -p $(dir $@); fi;)
		@$(ALICE_ROOT)/share/alibtool dependF "$(JETANELIBSDIR) $(JETANINC) $(DEPINC)  $<" > $@
$(MODDIRO)/%.d: $(MODDIRS)/%.F
ifndef ALIQUIET
		@echo "***** Making dependences for $< *****";
endif
		@(if [ ! -d '$(dir $@)' ]; then echo "***** Making directory $(dir $@) *****"; mkdir -p $(dir $@); fi;)
		$(MUTE)share/alibtool dependF "$(JETANELIBSDIR) $(JETANINC) $(DEPINC)  $<" > $@
$(MODDIRO)/%.d: $(MODDIRS)/%.c
ifndef ALIQUIET
		@echo "***** Making dependences for $< *****";
endif
		@(if [ ! -d '$(dir $@)' ]; then echo "***** Making directory $(dir $@) *****"; mkdir -p $(dir $@); fi;)
		@$(ALICE_ROOT)/share/alibtool depend "$(JETANDEFINE) $(JETANELIBSDIR) $(JETANINC) $(DEPINC) $<" > $@

.PRECIOUS: $(patsubst %.cxx,$(MODDIRO)/%.d,$(SRCS))
.PRECIOUS: $(patsubst %.c,$(MODDIRO)/%.d,$(CSRCS))
.PRECIOUS: $(patsubst %.F,$(MODDIRO)/%.d,$(patsubst %.f,$(MODDIRO)/%.d,$(FSRCS)))


#  Clean option for the whole module
clean:
ifndef ALIQUIET
		@echo "***** Cleaning JETAN *****"
endif
		$(MUTE)rm -rf $(MODDIRO)
		$(MUTE)rm -rf lib
		$(MUTE)rm -rf $(BINPATH)
