# $Id$
# Makefile template for the Alice HLT TPC display library

MODULE 				= AliHLTTPCDisplay

AM_CPPFLAGS			= @ALIROOT_CPPFLAGS@         \
				  @HOMER_CPPFLAGS@ \
				  -I$(top_srcdir)/BASE       \
				  -I$(srcdir)/../

# library definition
lib_LTLIBRARIES			=  libAliHLTTPCDisplay.la

# library sources
libAliHLTTPCDisplay_la_SOURCES	= AliHLTTPCDisplayMain.cxx \
				  AliHLTTPCDisplay.cxx  \
				  AliHLTTPCDisplay3D.cxx  \
				  AliHLTTPCDisplayCharge.cxx  \
				  AliHLTTPCDisplayFront.cxx  \
				  AliHLTTPCDisplayPad.cxx  \
				  AliHLTTPCDisplayPadRow.cxx  \
				  AliHLTTPCDisplayResiduals.cxx

# class header files, the link definition for the root dictionary
# will be created from the names of the header files
CLASS_HDRS			= AliHLTTPCDisplayMain.h  \
				  AliHLTTPCDisplay.h  \
				  AliHLTTPCDisplay3D.h  \
				  AliHLTTPCDisplayCharge.h  \
				  AliHLTTPCDisplayFront.h  \
				  AliHLTTPCDisplayPad.h  \
				  AliHLTTPCDisplayPadRow.h  \
				  AliHLTTPCDisplayResiduals.h



pkginclude_HEADERS		= $(CLASS_HDRS)

# version info for the library
libAliHLTTPCDisplay_la_LDFLAGS	= -L@ROOTLIBDIR@ \
				  @ROOTLIBS@ \
				  @ALIROOT_LDFLAGS@ \
				  @ALIROOT_LIBS@ \
				  @HOMER_LDFLAGS@ \
				  -version-info 1:0:0

# set the file name for the generated root dictionary
DICTCPP				=  AliHLTTPCDisplay-DICT.cxx
# add additional include files which are necessary for the compilation of the 
# dictionary files and which are not part of the CLASS_HDRS
DICTINCLUDE			= AliHLTTPCSpacePointData.h \
				  AliHLTTPCTrackArray.h
nodist_libAliHLTTPCDisplay_la_SOURCES  =  $(DICTCPP)

CLEANFILES			= 

include $(top_srcdir)/make.dict
