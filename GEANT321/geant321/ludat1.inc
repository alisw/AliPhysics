*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:20:35  cernlib
* Geant
*
*
#ifndef CERNLIB_GEANT321_LUDAT1_INC
#define CERNLIB_GEANT321_LUDAT1_INC
*
*
* ludat1.inc
*
#if !defined(CERNLIB_INCLUDE)
      COMMON/LUDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
      INTEGER MSTU,MSTJ
      REAL    PARU,PARJ
      SAVE /LUDAT1/
*
#endif
*
* ludat1.inc
*
#if defined(CERNLIB_INCLUDE)
*
#include "geant321/incbeg.inc"
     +    ludat1
#include "geant321/incend.inc"
#endif

#endif
