#$Id$
###############################
# Makefile for HLT code.      #
#                             #
# Author: Anders Vestbo,      #
#         Constantin Loizides #                    
###############################

LIBS = src hough comp misc
MAKEDIRS = $(LIBS) programs

default: so

all: 
	for i in $(MAKEDIRS) ; do cd $$i && $(MAKE) ; cd .. ; done

so:
	for i in $(LIBS) ; do cd $$i && $(MAKE) so ; cd .. ; done

static: 
	for i in $(LIBS) ; do cd $$i && $(MAKE) static ; cd .. ; done

libs: 
	for i in $(LIBS) ; do cd $$i && $(MAKE) libs ; cd .. ; done

clean:
	for i in $(MAKEDIRS) ; do cd $$i && $(MAKE) clean ; cd .. ; done

realclean: clean
	rm -rf $(ALIHLT_LIBDIR) 

print:
	@echo "ROOTSYS           = $(ROOTSYS)"
	@echo "ALICE_ROOT        = $(ALICE_ROOT)"
	@echo ""
	@echo "ALIHLT_BASEDIR    = $(ALIHLT_BASEDIR)"
	@echo "ALIHLT_TOPDIR     = $(ALIHLT_TOPDIR)"
	@echo "ALIHLT_LIBDIR     = $(ALIHLT_LIBDIR)"
	@echo "ALIHLT_MLUCDIR    = $(ALIHLT_MLUCDIR)"
	@echo "ALIHLT_USEPACKAGE = $(ALIHLT_USEPACKAGE)"
	@echo "ALIHLT_NOLOGGING  = $(ALIHLT_NOLOGGING)"
	@echo "ALIHLT_DOMC       = $(ALIHLT_DOMC)"
	@echo "ALIHLT_ALIDETECT  = $(ALIHLT_ALIDETECT)"

help:
	cat doc/README

.phony: help print clean realclean
