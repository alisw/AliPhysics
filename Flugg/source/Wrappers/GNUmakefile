# $Id$
# --------------------------------------------------------------
# GNUmakefile for materials library.  Gabriele Cosmo, 18/9/96.
# --------------------------------------------------------------

name := FggWrappers

ifndef FLUGGINSTALL
  FLUGGINSTALL = ../..
endif

include $(FLUGGINSTALL)/config/architecture.gmk

CPPFLAGS += -I$(G4BASE_G4)/global/management/include \
	    -I$(G4BASE_G4)/global/STLInterface \
            -I$(G4BASE_G4)/global/HEPRandom/include \
            -I$(G4BASE_G4)/global/HEPGeometry/include \
	    -I$(G4BASE_G4)/geometry/management/include \
	    -I$(G4BASE_G4)/geometry/volumes/include \
	    -I$(G4BASE_G4)/geometry/magneticfield/include \
	    -I$(G4BASE_G4)/materials/include \
	    -I$(G4BASE)/dummyG4/include \
	    -I$(G4BASE)/Wrappers/include
	    
include $(FLUGGINSTALL)/config/common.gmk

.PHONY: global
 
global: lib

