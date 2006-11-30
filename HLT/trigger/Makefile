#$Id$
##################################
# Makefile for HLT Trigger code. #
#                                #
# Author: Anders Vestbo,         #                    
#         Constantin Loizides    #
##################################

MODNAME = AliHLTTrigger

ifndef ALIHLT_TOPDIR
TOPDIR = $(shell pwd)/..
endif


include $(ALIHLT_TOPDIR)/Makefile.conf

SRCS	= AliD0toKpi.cxx AliD0Trigger.cxx 

include $(ALIHLT_TOPDIR)/Makefile.rules

INCLUDES += -I$(ALICE_ROOT)/ITS
