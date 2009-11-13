
AM_CPPFLAGS			= -I$(top_srcdir)/BASE \
				  -I$(srcdir)/.. \
				  -I$(top_srcdir)/MUON \
				  -I$(ALICE_ROOT)/MUON \
				  @ALIROOT_CPPFLAGS@ \
				  -I@ROOTINCDIR@

EXTRA_DIST			= 

check_PROGRAMS			= testTriggerDomain testGlobalTriggerComponent testMuonSpectroTrigger

testTriggerDomain_SOURCES 	= testTriggerDomain.C

testGlobalTriggerComponent_SOURCES	= testGlobalTriggerComponent.C

testMuonSpectroTrigger_SOURCES	= testMuonSpectroTrigger.C

# linker flags
LDADD_COMMON			= $(top_builddir)/BASE/libHLTbase.la \
				  $(top_builddir)/MUON/libAliHLTMUON.la \
				  $(top_builddir)/trigger/libAliHLTTrigger.la
LDFLAGS_COMMON			= -L@ROOTLIBDIR@ \
				  @ROOTLIBS@ \
				  @ALIROOT_LDFLAGS@ \
				  @ALIROOT_LIBS@

testTriggerDomain_LDADD 	= $(LDADD_COMMON)

testTriggerDomain_LDFLAGS	= $(LDFLAGS_COMMON)

testGlobalTriggerComponent_LDADD 	= $(LDADD_COMMON)

testGlobalTriggerComponent_LDFLAGS	= $(LDFLAGS_COMMON)

testMuonSpectroTrigger_LDADD 	= $(LDADD_COMMON)

testMuonSpectroTrigger_LDFLAGS	= $(LDFLAGS_COMMON)

TESTS				= $(check_PROGRAMS)

