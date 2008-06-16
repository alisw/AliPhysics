# $Id$
# Makefile template ALICE HLT libHLTrec library test programs

AM_CPPFLAGS			= -I$(top_srcdir)/BASE \
				  -I$(srcdir)/.. \
				  @ALIROOT_CPPFLAGS@ \
				  -I@ROOTINCDIR@

EXTRA_DIST			= 

check_PROGRAMS			= testAliHLTEsdManager \
				  testAliRawReaderFile

testAliHLTEsdManager_SOURCES 	= testAliHLTEsdManager.C

testAliRawReaderFile_SOURCES 	= testAliRawReaderFile.C


# linker flags
LDADD_COMMON		 	= $(top_builddir)/BASE/libHLTbase.la \
				  $(top_builddir)/rec/libHLTrec.la
LDFLAGS_COMMON			= -L@ROOTLIBDIR@ \
				  @ROOTLIBS@ \
				  @ALIROOT_LDFLAGS@ \
				  @ALIROOT_LIBS@

testAliHLTEsdManager_LDADD 	= $(LDADD_COMMON)

testAliHLTEsdManager_LDFLAGS	= $(LDFLAGS_COMMON)

testAliRawReaderFile_LDADD 	= $(LDADD_COMMON)

testAliRawReaderFile_LDFLAGS	= $(LDFLAGS_COMMON)

TESTS				= $(check_PROGRAMS)
