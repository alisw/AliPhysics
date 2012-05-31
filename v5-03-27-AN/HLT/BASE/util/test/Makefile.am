# $Id$
# Makefile template for Alice HLT libAliHLTUtil test programs

AM_CPPFLAGS			= -I$(top_srcdir)/BASE \
				  -I$(top_srcdir)/BASE/util \
				  @ALIROOT_CPPFLAGS@

EXTRA_DIST	= 

check_PROGRAMS	= testAliHLTFileWriter \
		  testAliHLTGlobalBarrelTrack

testAliHLTFileWriter_SOURCES = testAliHLTFileWriter.C
testAliHLTFileWriter_LDADD = -lHLTbase \
				    @ALIROOT_LDFLAGS@ \
				    @ALIROOT_LIBS@

testAliHLTGlobalBarrelTrack_SOURCES = testAliHLTGlobalBarrelTrack.C
testAliHLTGlobalBarrelTrack_LDADD = -lHLTbase \
				    -lAliHLTUtil \
				    @ALIROOT_LDFLAGS@ \
				    @ALIROOT_LIBS@

TESTS		= $(check_PROGRAMS)
