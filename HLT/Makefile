#$Id$
###########################
# Makefile for HLT code.  #
#                         #
# Author: Anders Vestbo   #                    
###########################

MAKEDIRS = src hough comp programs misc

all:
	for i in $(MAKEDIRS) ; do cd $$i && make ; cd .. ; done

clean:
	for i in $(MAKEDIRS) ; do cd $$i && make clean ; cd .. ; done