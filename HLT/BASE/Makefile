#$Id$
##########################################
# Makefile for AliRoot HLT framework.    #
#                                        #
# Author: Matthias Richter,              #                     
#         Timm Steinbeck                 #
#         Artur Szostack                 #
##########################################

MODNAME = HLTbase

ifndef ALIHLT_TOPDIR
ALIHLT_TOPDIR = $(shell pwd)/..
endif

include $(ALIHLT_TOPDIR)/Makefile.conf

SRCS	= AliHLTComponent.cxx AliHLTComponentHandler.cxx \
          AliHLTSystem.cxx AliHLT_C_Component_WrapperInterface.cxx \
	  AliHLTProcessor.cxx AliHLTConfiguration.cxx \
	  AliHLTLogging.cxx AliHLTDataBuffer.cxx \
	  AliHLTDataSource.cxx AliHLTDataSink.cxx

SRCS += $(EXTRA_SRCS)


include $(ALIHLT_TOPDIR)/Makefile.rules
