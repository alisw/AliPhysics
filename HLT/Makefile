#$Id$
###########################
# Makefile for HLT code.  #
#                         #
# Author: Anders Vestbo   #                    
###########################

LIBS = src hough comp misc
MAKEDIRS = $(LIBS) programs

all:
	for i in $(MAKEDIRS) ; do cd $$i && $(MAKE) ; cd .. ; done

libs: 
	for i in $(LIBS) ; do cd $$i && $(MAKE) ; cd .. ; done

clean:
	for i in $(MAKEDIRS) ; do cd $$i && $(MAKE) clean ; cd .. ; done
