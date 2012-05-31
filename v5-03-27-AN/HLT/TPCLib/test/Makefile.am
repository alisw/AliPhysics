# $Id$
# Makefile template Alice HLT TPC library test programs

AM_CPPFLAGS			= -I$(top_srcdir)/BASE \
				  -I$(top_srcdir)/RCU \
				  -I$(srcdir)/.. \
				  @ALIROOT_CPPFLAGS@ \
				  -I@ROOTINCDIR@

EXTRA_DIST			= 

check_PROGRAMS			= testAliHLTTPCDigitReaderDecoder \
				  testAliHLTTPCDigitReaderPacked \
				  testAliHLTTPCMapping \
				  testAliHLTTPCDefinitions \
			          testAliHLTTPCDataCheckerComponent

testAliHLTTPCDigitReaderDecoder_SOURCES = testAliHLTTPCDigitReaderDecoder.C
testAliHLTTPCDigitReaderPacked_SOURCES	= testAliHLTTPCDigitReaderPacked.C
testAliHLTTPCMapping_SOURCES		= testAliHLTTPCMapping.C
testAliHLTTPCDefinitions_SOURCES	= testAliHLTTPCDefinitions.C
testAliHLTTPCDataCheckerComponent_SOURCES	= testAliHLTTPCDataCheckerComponent.C


# linker flags
LDADD_COMMON			 	= $(top_builddir)/BASE/libHLTbase.la \
					  $(top_builddir)/TPCLib/libAliHLTTPC.la
LDFLAGS_COMMON				= -L@ROOTLIBDIR@ \
					  @ROOTLIBS@ \
					  @ALIROOT_LDFLAGS@ \
					  @ALIROOT_LIBS@

testAliHLTTPCDigitReaderDecoder_LDADD 	= $(LDADD_COMMON)
testAliHLTTPCDigitReaderPacked_LDADD 	= $(LDADD_COMMON)
testAliHLTTPCMapping_LDADD 		= $(LDADD_COMMON)
testAliHLTTPCDefinitions_LDADD		= $(LDADD_COMMON)
testAliHLTTPCDataCheckerComponent_LDADD	= $(LDADD_COMMON)
testAliHLTTPCDigitReaderDecoder_LDFLAGS	= $(LDFLAGS_COMMON)
testAliHLTTPCDigitReaderPacked_LDFLAGS	= $(LDFLAGS_COMMON)
testAliHLTTPCMapping_LDFLAGS		= $(LDFLAGS_COMMON)
testAliHLTTPCDefinitions_LDFLAGS	= $(LDFLAGS_COMMON)
testAliHLTTPCDataCheckerComponent_LDFLAGS	= $(LDFLAGS_COMMON)

# set back to all as sson as DigitReaderPacked is fixed
#TESTS				= $(check_PROGRAMS)
TESTS				= testAliHLTTPCDigitReaderDecoder \
				  testAliHLTTPCMapping \
				  testAliHLTTPCDefinitions \
				  testAliHLTTPCDataCheckerComponent
