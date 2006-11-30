#$Id$
################################
# Makefile for HLT Hough code. #
#                              #
# Author: Anders Vestbo,       #                    
#         Constantin Loizides  #
################################

MODNAME = AliHLTKalman

ifndef ALIHLT_TOPDIR
ALIHLT_TOPDIR = $(shell pwd)/..
endif

include $(ALIHLT_TOPDIR)/Makefile.conf

SRCS = AliHLTKalman.cxx AliHLTKalmanTrack.cxx 

ifeq ($(ALIHLT_USEPACKAGE),ROOT)
HDRS = $(SRCS:.cxx=.h) $(MODNAME)LinkDef.h
OBJS = $(patsubst %.cxx,$(OBJDIR)/%.o,$(SRCS)) $(DICTO)
endif

ifeq ($(ALIHLT_USEPACKAGE),ALIROOT)
HDRS = $(SRCS:.cxx=.h) $(MODNAME)LinkDef.h
OBJS = $(patsubst %.cxx,$(OBJDIR)/%.o,$(SRCS)) $(DICTO)
endif

include $(ALIHLT_TOPDIR)/Makefile.rules
