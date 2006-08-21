# $Id$
# Makefile template for the Alice HLT sample library

EXTRA_DIST			= tutorial.c

MODULE 				= AliHLTSample

AM_CPPFLAGS			= -I$(top_srcdir)/BASE


# library definition
lib_LTLIBRARIES			=  libAliHLTSample.la

# library sources
libAliHLTSample_la_SOURCES	=  AliHLTSampleComponent1.cxx \
				   AliHLTSampleComponent2.cxx \
				   AliHLTDummyComponent.cxx

# class header files, the link definition for the root dictionary
# will be created from the names of the header files
CLASS_HDRS			=  AliHLTSampleComponent1.h \
				   AliHLTSampleComponent2.h \
				   AliHLTDummyComponent.h

noinst_HEADERS			=  $(CLASS_HDRS)

# version info for the library
libAliHLTSample_la_LDFLAGS	=  -version-info 1:0:0

# set the file name for the generated root dictionary
DICTCPP				=  AliHLTSample-DICT.cxx
nodist_libAliHLTSample_la_SOURCES    =  $(DICTCPP)

CLEANFILES			=

include $(top_srcdir)/make.dict
