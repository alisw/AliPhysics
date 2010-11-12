
AM_CPPFLAGS			= -I$(top_srcdir)/BASE \
				  -I$(srcdir)/.. \
				  -I$(top_srcdir)/MUON \
				  -I$(ALICE_ROOT)/MUON \
				  @ALIROOT_CPPFLAGS@ \
				  -I@ROOTINCDIR@

EXTRA_DIST			= 

check_PROGRAMS			= testTriggerDomain testGlobalTriggerComponent testMuonSpectroTrigger testAliHLTTriggerCounters testTriggerCounterComponent

testTriggerDomain_SOURCES 	= testTriggerDomain.C

testGlobalTriggerComponent_SOURCES	= testGlobalTriggerComponent.C

testMuonSpectroTrigger_SOURCES	= testMuonSpectroTrigger.C

testAliHLTTriggerCounters_SOURCES = testAliHLTTriggerCounters.C

testTriggerCounterComponent_SOURCES	= testTriggerCounterComponent.C

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

testAliHLTTriggerCounters_LDADD   = $(LDADD_COMMON)
testAliHLTTriggerCounters_LDFLAGS = $(LDFLAGS_COMMON)

testTriggerCounterComponent_LDADD 	= $(LDADD_COMMON)
testTriggerCounterComponent_LDFLAGS	= $(LDFLAGS_COMMON)

TESTS				= $(check_PROGRAMS)

