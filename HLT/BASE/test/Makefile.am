# $Id$
# Makefile template Alice HLT framework test programs

AM_CPPFLAGS			= -I$(top_srcdir)/BASE

EXTRA_DIST	= 

check_PROGRAMS	= dtOperators \
		  testAliHLT_C_Component_WrapperInterface

dtOperators_SOURCES = dtOperators.cxx
dtOperators_LDADD = $(top_builddir)/BASE/libHLTbase.la

testAliHLT_C_Component_WrapperInterface_SOURCES = testAliHLT_C_Component_WrapperInterface.C
testAliHLT_C_Component_WrapperInterface_LDADD = $(top_builddir)/BASE/libHLTbase.la

TESTS		= $(check_PROGRAMS)