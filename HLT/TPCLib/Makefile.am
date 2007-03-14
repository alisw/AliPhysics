# $Id$
# Makefile template for the Alice HLT TPC library

MODULE 				= AliHLTTPC

EXTRA_DIST			= AliHLTTPCLinkDef.h

AM_CPPFLAGS			= @ALIROOT_CPPFLAGS@         \
				  -I$(top_srcdir)/BASE       \
				  -I$(top_srcdir)/BASE/util

# generation of the TPC Pad mapping tables for the 
# AliHLTTPCDigitReaderRaw reader
if USE_TPC_MAPPING
TPC_MAPPING_NOINST_SRC		= mapping2array.cxx
TPC_MAPPING_INC			= mapping_array_out.inc
TPC_MAPPING_NOINST_PGM		= mapping2array
TPC_MAPPING_COMPONENTS		= AliHLTTPCDigitReaderRaw.cxx
endif

noinst_PROGRAMS			= $(TPC_MAPPING_NOINST_PGM)
if USE_TPC_MAPPING
mapping2array_SOURCES		= $(TPC_MAPPING_NOINST_SRC)
endif

# library definition
lib_LTLIBRARIES			=  libAliHLTTPC.la

# version info for the library
LIBRARY_VERSION			= '3:0:0'

# MODDIR is set by the AliRoot build system and denotes the topdir
# of the module, we must set it since the package definition libAliHLTTPC.pkg
# includes another common configuration file
MODDIR				= $(top_srcdir)
PKGDEF				= $(MODDIR)/libHLTbase.pkg
include $(top_srcdir)/libAliHLTTPC.pkg

# library sources
libAliHLTTPC_la_SOURCES		= $(MODULE_SRCS)

# library headers
pkginclude_HEADERS		= $(MODULE_HDRS)

# linker flags
libAliHLTTPC_la_LIBADD		= $(top_builddir)/BASE/util/libAliHLTUtil.la
libAliHLTTPC_la_LDFLAGS		= -L@ROOTLIBDIR@ \
				  @ROOTLIBS@ \
				  @ALIROOT_LDFLAGS@ \
				  @ALIROOT_LIBS@ \
				  -version-info  $(LIBRARY_VERSION)

# set the file name for the generated root dictionary
DICTCPP				=  AliHLTTPC-DICT.cxx
nodist_libAliHLTTPC_la_SOURCES  =  $(DICTCPP)

BUILT_SOURCES			= $(TPC_MAPPING_INC)
CLEANFILES			= $(BUILT_SOURCES)

include $(top_srcdir)/make.dict

SUBDIRS 			= . OnlineDisplay

$(TPC_MAPPING_INC): $(TPC_MAPPING_NOINST_PGM)
	./$<
