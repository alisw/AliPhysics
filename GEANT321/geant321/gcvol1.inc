*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:20:34  cernlib
* Geant
*
*
#ifndef CERNLIB_GEANT321_GCVOL1_INC
#define CERNLIB_GEANT321_GCVOL1_INC
*
*
* gcvol1.inc
*
#if !defined(CERNLIB_INCLUDE)
      COMMON/GCVOL1/NLEVL1,NAMES1(15),NUMBR1(15),LVOLU1(15)
C
#if defined(CERNLIB_TYPE)
#include "geant321/gtvol1.inc"
#endif
#endif
*
* gcvol1.inc
*
#if defined(CERNLIB_INCLUDE)
*
#include "geant321/incbeg.inc"
     +   gcvol1
#include "geant321/incend.inc"
#endif

#endif
