#$Id$
################################
# Makefile for HLT Misc code.  #
#                              #
# Author: Anders Vestbo        #                    
#         Constantin Loizides  #
################################

MODNAME = AliL3Misc

ifndef ALIHLT_TOPDIR
ALIHLT_TOPDIR = $(shell pwd)/..
endif

include $(ALIHLT_TOPDIR)/Makefile.conf

SRCS	= AliL3TransBit.cxx AliL3AltroMemHandler.cxx \
          AliL3VHDLClusterFinder.cxx AliL3TPCMapping.cxx \
          AliL3FFloat.cxx AliL3Stopwatch.cxx AliL3DDLRawReader.cxx \
          AliL3DDLRawReaderFile.cxx AliL3DDLTPCRawStream.cxx 

ifneq ($(ARCH),macosx)
SRCS   += AliL3DDLDataFileHandler.cxx AliL3DataHandler.cxx 
endif

include $(ALIHLT_TOPDIR)/Makefile.rules
