*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:20:34  cernlib
* Geant
*
*
#ifndef CERNLIB_GEANT321_LUJETS_INC
#define CERNLIB_GEANT321_LUJETS_INC
*
*
* lujets.inc
*
#if !defined(CERNLIB_INCLUDE)
      COMMON/LUJETS/N,K(4000,5),P(4000,5),V(4000,5)
      INTEGER N,K
      REAL P,V
      SAVE /LUJETS/
*
#endif
*
* lujets.inc
*
#if defined(CERNLIB_INCLUDE)
*
#include "geant321/incbeg.inc"
     +    lujets
#include "geant321/incend.inc"
#endif

#endif
