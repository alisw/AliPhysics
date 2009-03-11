# $Id$
# Makefile template for the Alice HLT framework and components


if HAVE_ALIROOT
ALIROOT_DEP=sim \
	    rec \
	    shuttle \
	    pendolino
endif

if EN_HLT_SAMPLE
SAMPLE_DIR=SampleLib	
endif

if EN_HLT_TPC
TPC_DIR=TPCLib
endif

if EN_HLT_RCU
RCU_DIR=RCU
endif

if EN_HLT_PHOS
PHOS_DIR=PHOS
endif

if EN_HLT_TRD
TRD_DIR=TRD
endif

if EN_HLT_FMD
FMD_DIR=FMD
endif

if EN_HLT_MUON
MUON_DIR=MUON
endif

if EN_HLT_COMP
COMP_DIR=comp
endif

if EN_HLT_TRIGGER
TRIGGER_DIR=trigger
endif

if EN_HLT_ITS
ITS_DIR=ITS
endif

if EN_HLT_EMCAL
EMCAL_DIR=EMCAL
endif

if EN_HLT_GLOBAL
GLOBAL_DIR=global
endif

if EN_HLT_JET
JET_DIR=JET
endif


SUBDIRS 		= BASE \
			  $(ALIROOT_DEP) \
			  $(RCU_DIR) \
			  $(SAMPLE_DIR) \
			  $(TPC_DIR) \
			  $(PHOS_DIR) \
			  $(TRD_DIR) \
			  $(FMD_DIR) \
			  $(MUON_DIR) \
			  $(COMP_DIR) \
                          $(TRIGGER_DIR) \
                          $(ITS_DIR) \
                          $(EMCAL_DIR) \
                          $(GLOBAL_DIR) \
			  doc

#$(JET_DIR) \

EXTRA_DIST		= libHLTbase.pkg \
			  libHLTinterface.pkg \
			  libHLTsim.pkg \
			  libHLTrec.pkg \
			  libHLTshuttle.pkg \
			  libHLTpendolino.pkg \
			  libAliHLTHOMER.pkg \
			  libAliHLTUtil.pkg \
			  libAliHLTSample.pkg \
			  libAliHLTPHOS.pkg \
			  libAliHLTTPC.pkg \
			  libAliHLTTRD.pkg \
			  libAliHLTFMD.pkg \
			  libAliHLTMUON.pkg \
			  libAliHLTComp.pkg \
			  libAliHLTTrigger.pkg \
			  libAliHLTRCU.pkg \
			  libAliHLTITS.pkg \
			  libAliHLTEMCAL.pkg \
			  libAliHLTGlobal.pkg \
			  libAliHLTJet.pkg \
			  exa/sample-component1.C \
			  exa/monitoring.C \
			  hlt.conf \
			  .revision

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

svn-update:
	@cd $(top_srcdir) ;\
	revision=`svn info 2> /dev/null | grep "Revision:" | cut -d ' ' -f 2 | tr -d '\n'` ;\
	last=`svn info 2> /dev/null | grep "Last Changed Rev:" | cut -d ' ' -f 4 | cut -d '/' -f 4 | tr -d '\n'` ;\
	if test "x$$revision" != "x" && test $$last -gt $$revision; then \
		echo "The SVN copy is not up-to-date (revision $$revision - last update $$last), please 'svn update'"; \
		exit -1; \
	fi

.revision: svn-update
	@cd $(top_srcdir) ;\
	svn info 2> /dev/null | grep "Revision:" | cut -d ' ' -f 2 | tr -d '\n' > $@
#
# EOF
#
