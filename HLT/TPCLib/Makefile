#$Id$
##########################################
# Makefile for AliRoot HLT Sample        #
#                                        #
# Author: Matthias Richter,              #                     
##########################################

MODNAME = AliHLTTPC

ifndef ALIHLT_TOPDIR
ALIHLT_TOPDIR = $(shell pwd)/..
endif

USENEWIO=1

include $(ALIHLT_TOPDIR)/Makefile.conf

SRCS	= AliHLTTPCDefinitions.cxx \
          AliHLTTPCRawDataUnpackerComponent.cxx \
          AliHLTTPCClusterFinderComponent.cxx \
          AliHLTTPCVertexFinderComponent.cxx \
          AliHLTTPCSliceTrackerComponent.cxx \
          AliHLTTPCGlobalMergerComponent.cxx \
          AliRawReaderMemory.cxx

SRCS += $(EXTRA_SRCS)


include $(ALIHLT_TOPDIR)/Makefile.rules
