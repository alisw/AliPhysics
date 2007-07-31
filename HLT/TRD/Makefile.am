# $Id$
# Makefile template for the Alice HLT TRD library

MODULE 				= AliHLTTRD

EXTRA_DIST			= AliHLTTRDLinkDef.h

AM_CPPFLAGS			= @ALIROOT_CPPFLAGS@         \
				  -I$(top_srcdir)/BASE \
				  -I$(ALICE_ROOT)/TRD

# noinst_PROGRAMS			= $(TRD_MAPPING_NOINST_PGM)
# if USE_TRD_MAPPING
# mapping2array_SOURCES		= $(TRD_MAPPING_NOINST_SRC)
# endif

# library definition
lib_LTLIBRARIES			=  libAliHLTTRD.la

# library sources
libAliHLTTRD_la_SOURCES		= AliTRDclusterizerV1HLT.cxx \
				  AliTRDclusterizerV2HLT.cxx \
				  AliHLTTRDDefinitions.cxx \
				  AliHLTTRDClusterizerComponent.cxx \
				  AliHLTTRDClusterizerComponentV2.cxx \
				  AliHLTTRDTrackerComponent.cxx \
				  AliTRDtrackerHLT.cxx

# class header files, the link definition for the root dictionary
# will be created from the names of the header files
CLASS_HDRS			= AliTRDclusterizerV1HLT.h \
				  AliTRDclusterizerV2HLT.h \
				  AliHLTTRDDefinitions.h \
				  AliHLTTRDClusterizerComponent.h \
				  AliHLTTRDClusterizerComponentV2.h \
				  AliHLTTRDTrackerComponent.h \
				  AliTRDtrackerHLT.h

pkginclude_HEADERS		= $(CLASS_HDRS) 

# version info for the library
libAliHLTTRD_la_LDFLAGS		= -L@ROOTLIBDIR@ \
				  @ROOTLIBS@ \
				  @ALIROOT_LDFLAGS@ \
				  @ALIROOT_LIBS@ \
				  -version-info 1:0:0

# set the file name for the generated root dictionary
DICTCPP				=  AliHLTTRD-DICT.cxx
nodist_libAliHLTTRD_la_SOURCES  =  $(DICTCPP)

#BUILT_SOURCES			= $(TRD_MAPPING_INC)
#CLEANFILES			= $(BUILT_SOURCES)
CLEANFILES			= 
include $(top_srcdir)/make.dict

SUBDIRS 			= . #OnlineDisplay

#$(TRD_MAPPING_INC): $(TRD_MAPPING_NOINST_PGM)
#	./$<
