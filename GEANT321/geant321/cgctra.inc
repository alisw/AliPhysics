*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:20:42  cernlib
* Geant
*
*
#ifndef CERNLIB_GEANT321_CGCTRA_INC
#define CERNLIB_GEANT321_CGCTRA_INC
*
*
* cgctra.inc
*
#if !defined(CERNLIB_INCLUDE)
      PARAMETER (NTMAX = 3)
      COMMON /CGCTRA/ TSCRN(4,3,NTMAX),TBACK(4,3,NTMAX)
#endif
*
* cgctra.inc
*
#if defined(CERNLIB_INCLUDE)
#include "geant321/incbeg.inc"
       +   CGCTRA.
#include "geant321/incend.inc"
*
#endif

#endif
