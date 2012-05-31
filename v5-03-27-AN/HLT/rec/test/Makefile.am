# $Id$
# Makefile template ALICE HLT libHLTrec library test programs

AM_CPPFLAGS			= -I$(top_srcdir)/BASE \
				  -I$(srcdir)/.. \
				  @ALIROOT_CPPFLAGS@ \
				  -I@ROOTINCDIR@

EXTRA_DIST			= 

check_PROGRAMS			= testAliHLTEsdManager \
				  testAliRawReader

testAliHLTEsdManager_SOURCES 	= testAliHLTEsdManager.C

testAliRawReader_SOURCES 	= testAliRawReader.C


# linker flags
LDADD_COMMON		 	= -lHLTbase \
				  -lHLTrec
LDFLAGS_COMMON			= -L@ROOTLIBDIR@ \
				  @ROOTLIBS@ \
				  @ALIROOT_LDFLAGS@ \
				  @ALIROOT_LIBS@

testAliHLTEsdManager_LDADD 	= $(LDADD_COMMON)

testAliHLTEsdManager_LDFLAGS	= $(LDFLAGS_COMMON)

testAliRawReader_LDADD 		= $(LDADD_COMMON)

testAliRawReader_LDFLAGS	= $(LDFLAGS_COMMON)

TESTS				= $(check_PROGRAMS)
