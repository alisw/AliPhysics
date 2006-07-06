# $Id$
# Makefile template for the Alice HLT framework

MODULE				= HLTbase

AM_CPPFLAGS			= -I$(top_srcdir)/src      

# library definition
lib_LTLIBRARIES			=  libHLTbase.la

# library sources
libHLTbase_la_SOURCES		=  AliHLTComponent.cxx \
				   AliHLTComponentHandler.cxx \
				   AliHLTSystem.cxx \
				   AliHLT_C_Component_WrapperInterface.cxx \
				   AliHLTProcessor.cxx \
				   AliHLTConfiguration.cxx \
				   AliHLTLogging.cxx \
				   AliHLTDataBuffer.cxx \
				   AliHLTDataSource.cxx \
				   AliHLTDataSink.cxx

# class header files, the link definition for the root dictionary
# will be created from the names of the header files
CLASS_HDRS			=  AliHLTComponent.h \
				   AliHLTComponentHandler.h \
				   AliHLTSystem.h \
				   AliHLTProcessor.h \
				   AliHLTConfiguration.h \
				   AliHLTConfigurationHandler.h \
				   AliHLTTask.h \
				   AliHLTLogging.h \
				   AliHLTDataBuffer.h \
				   AliHLTDataSource.h \
				   AliHLTDataSink.h
pkginclude_HEADERS		=  $(CLASS_HDRS) \
				   AliHLTDataTypes.h \
				   AliHLT_C_Component_WrapperInterface.h \
				   AliHLTDefinitions.h


# version info for the library
libHLTbase_la_LDFLAGS		=  -version-info 1:0:0

# set the file name for the generated root dictionary
DICTCPP				=  HLTbase-DICT.cxx
nodist_libHLTbase_la_SOURCES    =  $(DICTCPP)

CLEANFILES			=

include $(top_srcdir)/make.dict
