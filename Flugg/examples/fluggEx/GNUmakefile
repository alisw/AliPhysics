# $Id$
# --------------------------------------------------------------
# GNUmakefile for examples module.  Gabriele Cosmo, 06/04/98.
# --------------------------------------------------------------

ifndef FLUGGINSTALL
  FLUGGINSTALL = ../..
endif
 
include $(FLUGGINSTALL)/config/architecture.gmk

SUBDIRS = AlAuAl BiasEx MagSphereNotRep T36flugg 

.PHONY : all clean clean_libs

all:
	@for dir in $(SUBDIRS); do (cd $$dir; $(MAKE)); done

clean:
	@for dir in $(SUBDIRS); do (cd $$dir; $(MAKE) clean); done

clean_libs:
	@for dir in $(SUBDIRS); do (cd $$dir; $(MAKE) clean_libs); done
