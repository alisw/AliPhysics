*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:20:42  cernlib
* Geant
*
*
#ifndef CERNLIB_GEANT321_CGCFAC_INC
#define CERNLIB_GEANT321_CGCFAC_INC
*
*
* cgcfac.inc
*
#if !defined(CERNLIB_INCLUDE)
      PARAMETER (LABC=1500)
      COMMON /CGCFAC/ F1(LABC),F2(LABC),F3(LABC),RFMIN(3),RFMAX(3),
     +                AABCD(4),DRFACE(5)
#endif
*
* cgcfac.inc
*
#if defined(CERNLIB_INCLUDE)
#include "geant321/incbeg.inc"
       +   CGCFAC.
#include "geant321/incend.inc"
*
#endif

#endif
