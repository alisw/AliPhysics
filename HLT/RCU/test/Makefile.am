# $Id$
# Makefile template Alice HLT RCU library test programs

AM_CPPFLAGS			= -I$(top_srcdir)/BASE \
				  -I$(srcdir)/.. \
				  @ALIROOT_CPPFLAGS@ \
				  -I@ROOTINCDIR@

EXTRA_DIST			= 

check_PROGRAMS			= testAliHLTAltroEncoder

testAliHLTAltroEncoder_SOURCES 	= testAliHLTAltroEncoder.C
testAliHLTAltroEncoder_LDADD 	= $(top_builddir)/BASE/libHLTbase.la \
				  $(top_builddir)/RCU/libAliHLTRCU.la

# linker flags
testAliHLTAltroEncoder_LDFLAGS	= -L@ROOTLIBDIR@ \
				  @ROOTLIBS@ \
				  @ALIROOT_LDFLAGS@ \
				  @ALIROOT_LIBS@

TESTS				= $(check_PROGRAMS)