*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:20:45  cernlib
* Geant
*
*
#ifndef CERNLIB_GEANT321_MMICAP_INC
#define CERNLIB_GEANT321_MMICAP_INC
*
*
* mmicap.inc
*
C
#include "geant321/mmicab.inc"
#include "geant321/gcbank.inc"
      DIMENSION LD(1),D(1)
      EQUIVALENCE (D(1),Q(1)),(LD(1),IQ(1))
      CHARACTER*24 DATSTR
      CHARACTER*80 COMMEN
      COMMON/ MICMAT / DATSTR,COMMEN,MATIDS(100,20,2)
C

#endif
