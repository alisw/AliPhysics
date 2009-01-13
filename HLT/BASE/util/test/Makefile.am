# $Id$
# Makefile template for Alice HLT libAliHLTUtil test programs

AM_CPPFLAGS			= -I$(top_srcdir)/BASE \
				  -I$(top_srcdir)/BASE/util

EXTRA_DIST	= 

check_PROGRAMS	= testAliHLTFileWriter

testAliHLTFileWriter_SOURCES = testAliHLTFileWriter.C
testAliHLTFileWriter_LDADD = $(top_builddir)/BASE/libHLTbase.la

TESTS		= $(check_PROGRAMS)
