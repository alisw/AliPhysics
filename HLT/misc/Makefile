#$Id$
################################
# Makefile for HLT Misc code.  #
#                              #
# Author: Anders Vestbo        #                    
#         Constantin Loizides  #
################################

MODNAME = AliHLTMisc

ifndef ALIHLT_TOPDIR
ALIHLT_TOPDIR = $(shell pwd)/..
endif

include $(ALIHLT_TOPDIR)/Makefile.conf

SRCS	= AliHLTTransBit.cxx AliHLTAltroMemHandler.cxx \
          AliHLTVHDLClusterFinder.cxx AliHLTTPCMapping.cxx \
          AliHLTFFloat.cxx AliHLTStopwatch.cxx AliHLTDDLRawReader.cxx \
          AliHLTDDLRawReaderFile.cxx AliHLTDDLTPCRawStream.cxx 

ifneq ($(ARCH),macosx)
## for macosx put into src as symbolic link
SRCS   += AliHLTDDLDataFileHandler.cxx AliHLTDataHandler.cxx 
endif

include $(ALIHLT_TOPDIR)/Makefile.rules
