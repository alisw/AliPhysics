# Sara Vanini 3/III/2000
# -----------------------------------------------------------
# GNUmakefile for dummyG4 classes library.  
# -----------------------------------------------------------

name := G4dummy

ifndef FLUGGINSTALL
  FLUGGINSTALL = ../..
endif

include $(FLUGGINSTALL)/config/architecture.gmk

CPPFLAGS += \
              -I$(G4BASE_G4)/global/management/include \
	      -I$(G4BASE_G4)/global/STLInterface \
              -I$(G4BASE_G4)/global/HEPRandom/include \
              -I$(G4BASE_G4)/global/HEPGeometry/include \
              -I$(G4BASE_G4)/global/HEPNumerics/include \
              -I$(G4BASE_G4)/geometry/management/include 

include $(FLUGGINSTALL)/config/common.gmk

.PHONY: global
 
global: lib

