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

SRCS	= AliL3TransBit.cxx AliL3AltroMemHandler.cxx AliL3DataHandler.cxx \
          AliL3VHDLClusterFinder.cxx AliL3TPCMapping.cxx \
          AliL3FFloat.cxx AliL3Stopwatch.cxx

include $(ALIHLT_TOPDIR)/Makefile.rules
