
AM_CPPFLAGS			= -I$(top_srcdir)/BASE \
				  -I$(srcdir)/.. \
				  @ALIROOT_CPPFLAGS@ \
				  -I@ROOTINCDIR@

EXTRA_DIST			= 

check_PROGRAMS			= testTriggerDomain

testTriggerDomain_SOURCES 	= testTriggerDomain.C


# linker flags
LDADD_COMMON			= $(top_builddir)/BASE/libHLTbase.la \
				  $(top_builddir)/trigger/libAliHLTTrigger.la
LDFLAGS_COMMON			= -L@ROOTLIBDIR@ \
				  @ROOTLIBS@ \
				  @ALIROOT_LDFLAGS@ \
				  @ALIROOT_LIBS@

testTriggerDomain_LDADD 	= $(LDADD_COMMON)

testTriggerDomain_LDFLAGS	= $(LDFLAGS_COMMON)

TESTS				= $(check_PROGRAMS)

