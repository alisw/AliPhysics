*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:20:35  cernlib
* Geant
*
*
#ifndef CERNLIB_GEANT321_GCPHRD_INC
#define CERNLIB_GEANT321_GCPHRD_INC
*
*
* gcphrd.inc
*
#if !defined(CERNLIB_INCLUDE)
C  GRATE - radiative modes' rates
      INTEGER KSHLS,L1SHLS,L2SHLS,L3SHLS,ISHLS,ISHLUS,ISHLTR
      REAL GPHRAT
      PARAMETER (KSHLS=6)
      PARAMETER (L1SHLS=8)
      PARAMETER (L2SHLS=7)
      PARAMETER (L3SHLS=8)
      PARAMETER (ISHLS=29)
      COMMON / GCPHRD / GPHRAT(ISHLS,MAXELZ),ISHLUS(24,4),ISHLTR(ISHLS)
C
#endif
*
* gcphrd.inc
*
#if defined(CERNLIB_INCLUDE)
*
#include "geant321/incbeg.inc"
     +   gcphrd
#include "geant321/incend.inc"
#endif

#endif
