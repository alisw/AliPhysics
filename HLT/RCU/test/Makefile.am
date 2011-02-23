# $Id$
# Makefile template Alice HLT RCU library test programs

AM_CPPFLAGS			= -I$(top_srcdir)/BASE \
				  -I$(srcdir)/.. \
				  @ALIROOT_CPPFLAGS@ \
				  -I@ROOTINCDIR@

EXTRA_DIST			= 

check_PROGRAMS			= testAliHLTAltroEncoder \
				  testAliHLTAltroGenerator

testAliHLTAltroEncoder_SOURCES 	= testAliHLTAltroEncoder.C
testAliHLTAltroGenerator_SOURCES= testAliHLTAltroGenerator.C


# linker flags
LDADD_COMMON		 	= -lHLTbase \
				  -lAliHLTRCU
LDFLAGS_COMMON			= -L@ROOTLIBDIR@ \
				  @ROOTLIBS@ \
				  @ALIROOT_LDFLAGS@ \
				  @ALIROOT_LIBS@

testAliHLTAltroEncoder_LDADD 	= $(LDADD_COMMON)
testAliHLTAltroGenerator_LDADD 	= $(LDADD_COMMON)

testAliHLTAltroEncoder_LDFLAGS	= $(LDFLAGS_COMMON)
testAliHLTAltroGenerator_LDFLAGS= $(LDFLAGS_COMMON)

TESTS				= $(check_PROGRAMS)
