*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:20:35  cernlib
* Geant
*
*
#ifndef CERNLIB_GEANT321_GC10EV_INC
#define CERNLIB_GEANT321_GC10EV_INC
*
*
* gc10ev.inc
*
#if !defined(CERNLIB_INCLUDE)
      REAL G10EV,TENEV
      PARAMETER (G10EV=1.0E-8)
      PARAMETER (TENEV=1.E-2)
C
#endif
*
* gc10ev.inc
*
#if defined(CERNLIB_INCLUDE)
*
#include "geant321/incbeg.inc"
     +   gc10ev
#include "geant321/incend.inc"
#endif

#endif
