*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:20:35  cernlib
* Geant
*
*
#ifndef CERNLIB_GEANT321_GCLUND_INC
#define CERNLIB_GEANT321_GCLUND_INC
*
*
* gclund.inc
*
#if !defined(CERNLIB_INCLUDE)
      COMMON/GCLUND/IFLUND,ECLUND
      INTEGER IFLUND
      REAL ECLUND
C
#endif
*
* gclund.inc
*
#if defined(CERNLIB_INCLUDE)
*
#include "geant321/incbeg.inc"
     +   gclund
#include "geant321/incend.inc"
#endif

#endif
