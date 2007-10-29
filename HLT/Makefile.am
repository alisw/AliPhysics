# $Id$
# Makefile template for the Alice HLT framework and components


if EN_HLT_SAMPLE
SAMPLE_DIR=SampleLib	
endif

if EN_HLT_TPC
TPC_DIR=TPCLib	
endif

if EN_HLT_PHOS
PHOS_DIR=PHOS
endif

if EN_HLT_TRD
TRD_DIR=TRD
endif

if EN_HLT_MUON
MUON_DIR=MUON
endif

if EN_HLT_COMP
COMP_DIR=comp
endif

SUBDIRS 		= BASE \
			  sim \
			  rec \
			  $(SAMPLE_DIR) \
			  $(TPC_DIR) \
			  $(PHOS_DIR) \
			  $(TRD_DIR) \
			  $(MUON_DIR) \
			  $(COMP_DIR) \
			  doc

EXTRA_DIST		= libHLTbase.pkg \
			  libHLTinterface.pkg \
			  libHLTsim.pkg \
			  libHLTrec.pkg \
			  libAliHLTHOMER.pkg \
			  libAliHLTUtil.pkg \
			  libAliHLTSample.pkg \
			  libAliHLTPHOS.pkg \
			  libAliHLTTPC.pkg \
			  libAliHLTTRD.pkg \
			  libAliHLTMUON.pkg \
			  libAliHLTComp.pkg \
			  hlt.conf

DIST_SUBDIRS		= $(SUBDIRS)

# this is a special target to create a stand-alone package from the SampleLib
SampleLibPkg: Makefile $(top_srcdir)/SampleLib/*
	@tmpdir=/tmp/$(USER)/$@ ; \
	if test -d $@ ; then \
	  cp -ruv $(top_srcdir)/SampleLib/* $@; \
	  cp -uv $(top_srcdir)/libAliHLTSample.pkg $@; \
	  cp -uv $(top_srcdir)/acinclude.m4 $@; \
	else \
	  (test ! -d $$tmpdir || rm $$tmpdir); \
	  echo "creating package in tmp dir $$tmpdir" ; \
	  mkdir -p $$tmpdir && \
	  (cp -r $(top_srcdir)/SampleLib/* $$tmpdir ; \
	  cp -uv $(top_srcdir)/libAliHLTSample.pkg $$tmpdir ; \
	  cp -uv $(top_srcdir)/acinclude.m4 $$tmpdir; \
	  (echo 'running autoreconf -f -i' ;cd $$tmpdir ; autoreconf -f -i) ; \
	  echo moving $$tmpdir to `pwd`; mv $$tmpdir . ) \
	fi

clean-local: 
	(test ! -d SampleLibPkg || rm -rf SampleLibPkg)

#
# EOF
#
