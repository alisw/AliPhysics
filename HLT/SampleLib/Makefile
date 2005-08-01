#$Id$
##########################################
# Makefile for AliRoot HLT Sample        #
#                                        #
# Author: Matthias Richter,              #                     
##########################################

MODNAME = AliHLTSample

ifndef ALIHLT_TOPDIR
ALIHLT_TOPDIR = $(shell pwd)/..
endif

include $(ALIHLT_TOPDIR)/Makefile.conf

SRCS	= AliHLTSampleComponent1.cxx AliHLTSampleComponent2.cxx 

SRCS += $(EXTRA_SRCS)


include $(ALIHLT_TOPDIR)/Makefile.rules
