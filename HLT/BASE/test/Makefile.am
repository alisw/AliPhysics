# $Id$
# Makefile template Alice HLT framework test programs

AM_CPPFLAGS			= -I$(top_srcdir)/BASE

EXTRA_DIST	= 

check_PROGRAMS	= testAliHLTBlockDataCollection \
		  testAliHLTComponent \
		  testAliHLTCTPData \
		  testAliHLTDataBuffer \
		  testAliHLTReadoutList \
		  testAliHLTScalars \
		  dtOperators \
		  testDefaultDataTypes

# 		  testAliHLT_C_Component_WrapperInterface

testAliHLTBlockDataCollection_SOURCES = testAliHLTBlockDataCollection.cxx
testAliHLTBlockDataCollection_LDFLAGS = @ALIROOT_LDFLAGS@ \
				      -lHLTbase

dtOperators_SOURCES = dtOperators.cxx
dtOperators_LDFLAGS = 	@ALIROOT_LDFLAGS@ \
			-lHLTbase

# testAliHLT_C_Component_WrapperInterface_SOURCES = testAliHLT_C_Component_WrapperInterface.C
# testAliHLT_C_Component_WrapperInterface_LDFLAGS = @ALIROOT_LDFLAGS@ \
						  -lHLTbase

testDefaultDataTypes_SOURCES = testDefaultDataTypes.C
testDefaultDataTypes_LDFLAGS = @ALIROOT_LDFLAGS@ \
			-lHLTbase

testAliHLTComponent_SOURCES = testAliHLTComponent.C
testAliHLTComponent_LDFLAGS = @ALIROOT_LDFLAGS@ \
			    -lHLTbase

testAliHLTCTPData_SOURCES = testAliHLTCTPData.C
testAliHLTCTPData_LDFLAGS = @ALIROOT_LDFLAGS@ \
			  -lHLTbase \
			  -lHLTrec

testAliHLTDataBuffer_SOURCES = testAliHLTDataBuffer.C
testAliHLTDataBuffer_LDFLAGS = @ALIROOT_LDFLAGS@ \
			    -lHLTbase

testAliHLTReadoutList_SOURCES = testAliHLTReadoutList.C
testAliHLTReadoutList_LDFLAGS = -lHLTbase \
				-L@ROOTLIBDIR@ \
				@ROOTLIBS@ \
				@ALIROOT_LDFLAGS@ \
				@ALIROOT_LIBS@

testAliHLTScalars_SOURCES = testAliHLTScalars.C
testAliHLTScalars_LDFLAGS = -lHLTbase \
			  -L@ROOTLIBDIR@ \
			  @ROOTLIBS@ \
			  @ALIROOT_LDFLAGS@ \
			  @ALIROOT_LIBS@

testAliHLTEventDDLBackwardCompatibility_SOURCES = testAliHLTEventDDLBackwardCompatibility.C
testAliHLTEventDDLBackwardCompatibility_LDFLAGS = -lHLTbase \
				-L@ROOTLIBDIR@ \
				@ROOTLIBS@ \
				@ALIROOT_LDFLAGS@ \
				@ALIROOT_LIBS@

TESTS		= $(check_PROGRAMS)
