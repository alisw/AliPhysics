# $Id$
# Makefile template Alice HLT framework test programs

AM_CPPFLAGS			= -I$(top_srcdir)/BASE

EXTRA_DIST	= 

check_PROGRAMS	= testAliHLTBlockDataCollection \
		  testAliHLTComponent \
		  testAliHLTCTPData \
		  testAliHLTDataBuffer \
		  testAliHLTReadoutList \
		  testAliHLTEventDDLBackwardCompatibility \
		  dtOperators \
		  testDefaultDataTypes

# 		  testAliHLT_C_Component_WrapperInterface

testAliHLTBlockDataCollection_SOURCES = testAliHLTBlockDataCollection.cxx
testAliHLTBlockDataCollection_LDADD = $(top_builddir)/BASE/libHLTbase.la

dtOperators_SOURCES = dtOperators.cxx
dtOperators_LDADD = $(top_builddir)/BASE/libHLTbase.la

# testAliHLT_C_Component_WrapperInterface_SOURCES = testAliHLT_C_Component_WrapperInterface.C
# testAliHLT_C_Component_WrapperInterface_LDADD = $(top_builddir)/BASE/libHLTbase.la

testDefaultDataTypes_SOURCES = testDefaultDataTypes.C
testDefaultDataTypes_LDADD = $(top_builddir)/BASE/libHLTbase.la

testAliHLTComponent_SOURCES = testAliHLTComponent.C
testAliHLTComponent_LDADD = $(top_builddir)/BASE/libHLTbase.la

testAliHLTCTPData_SOURCES = testAliHLTCTPData.C
testAliHLTCTPData_LDADD = $(top_builddir)/BASE/libHLTbase.la \
			  $(top_builddir)/rec/libHLTrec.la

testAliHLTDataBuffer_SOURCES = testAliHLTDataBuffer.C
testAliHLTDataBuffer_LDADD = $(top_builddir)/BASE/libHLTbase.la

testAliHLTReadoutList_SOURCES = testAliHLTReadoutList.C
testAliHLTReadoutList_LDADD = $(top_builddir)/BASE/libHLTbase.la \
				-L@ROOTLIBDIR@ \
				@ROOTLIBS@ \
				@ALIROOT_LDFLAGS@ \
				@ALIROOT_LIBS@

testAliHLTEventDDLBackwardCompatibility_SOURCES = testAliHLTEventDDLBackwardCompatibility.C
testAliHLTEventDDLBackwardCompatibility_LDADD = $(top_builddir)/BASE/libHLTbase.la \
				-L@ROOTLIBDIR@ \
				@ROOTLIBS@ \
				@ALIROOT_LDFLAGS@ \
				@ALIROOT_LIBS@

TESTS		= $(check_PROGRAMS)
