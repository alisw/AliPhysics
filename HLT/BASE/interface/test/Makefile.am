# $Id$
# Makefile template Alice HLT external interface programs

AM_CPPFLAGS			= -I$(top_srcdir)/BASE \
				  -I$(top_srcdir)/BASE/interface \
				  -DBASE_PATH=\"$(top_builddir)/BASE/interface/test/\" \
				  -I$(ALICE_ROOT)/STEER

EXTRA_DIST			= 

check_PROGRAMS			= testAliHLTExternalInterface

check_LTLIBRARIES		= libAliHLTTest.la

testAliHLTExternalInterface_SOURCES 	= testAliHLTExternalInterface.C \
					  $(top_srcdir)/BASE/AliHLTDataTypes.cxx

# Matthias 2009-05-29 at some point I want to get rid of the
# linking to libHLTbase. It should be enough just to load the
# interface library dynamically
# currently, the symbols of the HLTbase library are not recognized
# in the subsequent library loads if libHLTbase is loaded dynamically.
# Maybe an issue with the RTLD_NOW used in the dlopen from the
# AliHLTComponentHandler
testAliHLTExternalInterface_LDADD = -lHLTbase
testAliHLTExternalInterface_LDFLAGS  = -L@ROOTLIBDIR@ \
				  	@ROOTLIBS@ \
				  	@ALIROOT_LDFLAGS@ \
				  	@ALIROOT_LIBS@

libAliHLTTest_la_SOURCES = AliHLTTestProcessor.cxx AliHLTTestAgent.cxx
libAliHLTTest_la_LIBADD	 =
# the -rpath is the important thing here to build the shared
# library. Otherwise autotools just consider the libraries
# with the check_ attribute just convenience libraries which
# are static
libAliHLTTest_la_LDFLAGS = -rpath /tmp

TESTS		= $(check_PROGRAMS)
