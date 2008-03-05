# $Id$
# Makefile template Alice HLT framework test programs

AM_CPPFLAGS			= -I$(top_srcdir)/BASE

EXTRA_DIST	= 

check_PROGRAMS	= dtOperators

dtOperators_SOURCES = dtOperators.cxx
dtOperators_LDADD = $(top_builddir)/BASE/libHLTbase.la

TESTS		= $(check_PROGRAMS)