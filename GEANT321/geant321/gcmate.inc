*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:20:32  cernlib
* Geant
*
*
#ifndef CERNLIB_GEANT321_GCMATE_INC
#define CERNLIB_GEANT321_GCMATE_INC
*
*
* gcmate.inc
*
#if !defined(CERNLIB_INCLUDE)
      COMMON/GCMATE/NMAT,NAMATE(5),A,Z,DENS,RADL,ABSL
C
      INTEGER NMAT,NAMATE
      REAL A,Z,DENS,RADL,ABSL
C
#endif
*
* gcmate.inc
*
#if defined(CERNLIB_INCLUDE)
*
#include "geant321/incbeg.inc"
     +   gcmate
#include "geant321/incend.inc"
#endif

#endif
