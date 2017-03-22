#
# Libraries Makefile. Some ideas from Geant4 Makefiles
#
#                  M. Kirsanov 07.04.2006
#                     Modified 18.11.2006
#                     26.03.2008 CLHEP dependency removed
#                  N. Lavesson 28.04.2009 clean/distclean separated
#                  M. Kirsanov 21.07.2009 Mac-OSX flags added

SHELL = /bin/sh

-include config.mk

# flags:
#
#FFLAGSSHARED = -fPIC
CFLAGSSHARED = -fPIC
CXXFLAGSSHARED = -fPIC


HEPMCERROR=
ifneq (x$(HEPMCLOCATION),x)
 ifeq ($(wildcard $(HEPMCLOCATION)/include/HepMC/*.h),)
  HEPMCERROR= HepMC interface: ERROR, no HepMC headers found in ${HEPMCLOCATION}/include/HepMC
 endif
endif

# Location of directories.
MYTMPDIR=tmp
TOPDIR=$(shell \pwd)
INCDIR=include
SRCDIR=src
LIBDIR=lib
LIBDIRARCH=lib/archive
BINDIR=bin

# Location of libraries to be built.
ifeq ($(SHAREDLIBS),yes)
  targets=$(LIBDIRARCH)/libpythia8.a
  targets+=$(LIBDIR)/libpythia8.$(SHAREDSUFFIX)
  targets+=$(LIBDIRARCH)/liblhapdfdummy.a
  targets+=$(LIBDIR)/liblhapdfdummy.$(SHAREDSUFFIX)
else
  targets=$(LIBDIRARCH)/libpythia8.a
  targets+=$(LIBDIRARCH)/liblhapdfdummy.a
endif

ifneq (x$(HEPMCLOCATION),x)
 targets+=$(LIBDIRARCH)/libhepmcinterface.a
 ifeq ($(SHAREDLIBS),yes)
  targets+=$(LIBDIR)/libhepmcinterface.$(SHAREDSUFFIX)
 endif
endif


all: $(targets) config.mk

config.mk: ./configure
	./configure

# Main part: build Pythia8 library. 

$(MYTMPDIR)/%.o : $(SRCDIR)/%.cc
	@mkdir -p $(MYTMPDIR)
	$(CXX) $(CXXFLAGS) $(CXXFLAGSSHARED) -c -I$(INCDIR) $< -o $@

$(MYTMPDIR)/archive/%.o : $(SRCDIR)/%.cc
	@mkdir -p $(MYTMPDIR)/archive
	$(CXX) $(CXXFLAGS) -c -I$(INCDIR) $< -o $@

$(MYTMPDIR)/%.o : lhapdfdummy/%.cc
	@mkdir -p $(MYTMPDIR)
	$(CXX) $(CXXFLAGS) $(CXXFLAGSSHARED) -c -I$(INCDIR) $< -o $@

$(MYTMPDIR)/archive/%.o : lhapdfdummy/%.cc
	@mkdir -p $(MYTMPDIR)/archive
	$(CXX) $(CXXFLAGS) -c -I$(INCDIR) $< -o $@

# Creating the dependency files *.d
# The compiler with option -M is used to build the dependency strings. They
# are further edited with sed (stream editor). The first sed command adds the
# dependency for the *.d files themselves, the second one is needed because
# object files are put in the directory different from src. The last line
# removes empty *.d files produced in case of error.

ifeq ($(SHAREDLIBS),yes)
  $(MYTMPDIR)/%.d : $(SRCDIR)/%.cc
	@echo Making dependency for file $<; \
	mkdir -p $(MYTMPDIR); \
	$(CC) -M -I$(INCDIR) $< | \
	sed 's,\($*\)\.o[ :]*,\1.o $@ : ,g' | \
	sed 's/$*\.o/$(MYTMPDIR)\/$*.o/' > $@; \
	[ -s $@ ] || rm -f $@
endif

$(MYTMPDIR)/archive/%.d : $(SRCDIR)/%.cc
	@echo Making dependency for file $<; \
	mkdir -p $(MYTMPDIR)/archive; \
	$(CC) -M -I$(INCDIR) $< | \
	sed 's,\($*\)\.o[ :]*,\1.o $@ : ,g' | \
	sed 's/$*\.o/$(MYTMPDIR)\/archive\/$*.o/' > $@; \
	[ -s $@ ] || rm -f $@

objects := $(patsubst $(SRCDIR)/%.cc,$(MYTMPDIR)/%.o,$(wildcard $(SRCDIR)/*.cc))
objectsarch := $(patsubst $(SRCDIR)/%.cc,$(MYTMPDIR)/archive/%.o,$(wildcard $(SRCDIR)/*.cc))

$(LIBDIR)/libpythia8.$(SHAREDSUFFIX): $(objects)
	@mkdir -p $(LIBDIR)
	$(CXX) $(LDFLAGSSHARED) -o $@ $(objects) $(LDFLAGLIBNAME),$(notdir $@)

$(LIBDIRARCH)/libpythia8.a: $(objectsarch)
	@mkdir -p $(LIBDIRARCH)
	ar cru $@ $(objectsarch)

objdum := $(patsubst lhapdfdummy/%.cc,$(MYTMPDIR)/%.o,$(wildcard lhapdfdummy/*.cc))
objdumarch := $(patsubst lhapdfdummy/%.cc,$(MYTMPDIR)/archive/%.o,$(wildcard lhapdfdummy/*.cc))

$(LIBDIR)/liblhapdfdummy.$(SHAREDSUFFIX): $(objdum)
	@mkdir -p $(LIBDIR)
	$(CXX) $(LDFLAGSSHARED) -o $@ $(objdum) $(LDFLAGLIBNAME),$(notdir $@)

$(LIBDIRARCH)/liblhapdfdummy.a: $(objdumarch)
	@mkdir -p $(LIBDIRARCH)
	ar cru $@ $(objdumarch)

deps := $(patsubst $(SRCDIR)/%.cc,$(MYTMPDIR)/%.d,$(wildcard $(SRCDIR)/*.cc))
depsarch := $(patsubst $(SRCDIR)/%.cc,$(MYTMPDIR)/archive/%.d,$(wildcard $(SRCDIR)/*.cc))


# The "if" below is needed in order to avoid producing the dependency files
# when you want to just clean

ifeq (,$(findstring clean, $(MAKECMDGOALS)))
-include $(deps)
-include $(depsarch)
endif

# Build HepMC interface part if HepMC location is set.

ifneq (x$(HEPMCLOCATION),x)
 HEPMCINCLUDE=-I$(HEPMCLOCATION)/include

 ifeq (x$(HEPMCERROR),x)

   $(MYTMPDIR)/%.o : hepmcinterface/%.cc config.mk
	@mkdir -p $(MYTMPDIR)
	$(CXX) $(CXXFLAGS) $(CXXFLAGSSHARED) $(HEPMCVFLAG) -c -I$(INCDIR) $(HEPMCINCLUDE) $< -o $@

   $(MYTMPDIR)/archive/%.o : hepmcinterface/%.cc config.mk
	@mkdir -p $(MYTMPDIR)/archive
	$(CXX) $(CXXFLAGS) $(HEPMCVFLAG) -c -I$(INCDIR) $(HEPMCINCLUDE) $< -o $@

   $(MYTMPDIR)/%.d : hepmcinterface/%.cc
	@echo Making dependency for file $<; \
	mkdir -p $(MYTMPDIR); \
	$(CC) -M -I$(INCDIR) $(HEPMCINCLUDE) $< | \
	sed 's,\($*\)\.o[ :]*,\1.o $@ : ,g' | \
	sed 's/$*.o/$(MYTMPDIR)\/$*.o/' > $@; \
	[ -s $@ ] || rm -f $@

   $(MYTMPDIR)/archive/%.d : hepmcinterface/%.cc
	@echo Making dependency for file $<; \
	mkdir -p $(MYTMPDIR)/archive; \
	$(CC) -M -I$(INCDIR) $(HEPMCINCLUDE) $< | \
	sed 's,\($*\)\.o[ :]*,\1.o $@ : ,g' | \
	sed 's/$*.o/$(MYTMPDIR)\/archive\/$*.o/' > $@; \
	[ -s $@ ] || rm -f $@

   objectsI := $(patsubst hepmcinterface/%.cc,$(MYTMPDIR)/%.o,$(wildcard hepmcinterface/*.cc))
   objectsIarch := $(patsubst hepmcinterface/%.cc,$(MYTMPDIR)/archive/%.o,$(wildcard hepmcinterface/*.cc))

   $(LIBDIR)/libhepmcinterface.$(SHAREDSUFFIX) : $(objectsI)
	@mkdir -p $(LIBDIR)
	$(CXX) $(LDFLAGSSHARED) $(objectsI) -o $@ $(LDFLAGLIBNAME),$(notdir $@)

   $(LIBDIRARCH)/libhepmcinterface.a : $(objectsIarch)
	@mkdir -p $(LIBDIRARCH)
	ar cru $(LIBDIRARCH)/libhepmcinterface.a $(objectsIarch)

   depsI := $(patsubst hepmcinterface/%.cc,$(MYTMPDIR)/%.d,$(wildcard hepmcinterface/*.cc))
   depsIarch := $(patsubst hepmcinterface/%.cc,$(MYTMPDIR)/archive/%.d,$(wildcard hepmcinterface/*.cc))

   ifeq (,$(findstring clean, $(MAKECMDGOALS)))
   -include $(depsI)
   -include $(depsIarch)
   endif

 else

   $(LIBDIRARCH)/libhepmcinterface.a $(LIBDIR)/libhepmcinterface.$(SHAREDSUFFIX) :
	@echo $(HEPMCERROR)



 endif

endif


# Install targets:

ifneq (x$(INSTALLDIR),x.)
 install: all
	mkdir -p $(INSTALLDIR)
	mkdir -p $(DATADIR)
	make installit
 installit: installmain installdata
 installmain:
	cp -r include $(INSTALLDIR)/.
	cp -r lib $(INSTALLDIR)/.
	cp -p config.mk $(INSTALLDIR)/.

 ifneq ($(DATADIR),$(INSTALLDIR))
  installdata:
	rm -rf $(DATADIR)/xmldoc
	rm -f $(INSTALLDIR)/xmldoc
	cp -r *doc $(DATADIR)/.
	ln -fs $(DATADIR)/xmldoc $(INSTALLDIR)/xmldoc
 else
  installdata:
	rm -rf $(INSTALLDIR)/xmldoc
	cp -r xmldoc $(INSTALLDIR)/.
 endif
else
 install: all
	@echo "everything is already installed"
endif


# Clean up: remove (almost?) everything that cannot be recreated.

.PHONY: clean distclean

clean:
	rm -rf $(MYTMPDIR)
	rm -rf $(LIBDIR)
	rm -rf $(BINDIR)
	cd examples; rm -rf bin; rm -f *.exe; cd -
	cd rootexamples; rm -f *.exe; cd -

distclean: clean
	rm -f config.mk
	rm -f *~; rm -f \#*;
	cd $(SRCDIR); rm -f *~; rm -f \#*; cd -
	cd $(INCDIR); rm -f *~; rm -f \#*; cd -
	cd xmldoc; rm -f *~; rm -f \#*; cd -
	cd htmldoc; rm -f *~; rm -f \#*; cd -
	cd phpdoc; rm -f *~; rm -f \#*; cd -
	cd hepmcinterface; rm -f *~; rm -f \#*; cd -
	cd lhapdfdummy; rm -f *~; rm -f \#*; cd -
	cd examples; rm -f *~; rm -f \#*; rm -f core*; rm -f config.*; cd -
	cd rootexamples; rm -f *~; rm -f \#*; rm -f core*; rm -f config.*; cd -

