*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:20:35  cernlib
* Geant
*
*
#ifndef CERNLIB_GEANT321_GCPHPR_INC
#define CERNLIB_GEANT321_GCPHPR_INC
*
*
* gcphpr.inc
*
#if !defined(CERNLIB_INCLUDE)
C  Probability of radiative decay mode.
      REAL GFLUPR
      COMMON /GCPHPR/ GFLUPR(4,MAXELZ)
C
#endif
*
* gcphpr.inc
*
#if defined(CERNLIB_INCLUDE)
*
#include "geant321/incbeg.inc"
     +   gcphpr
#include "geant321/incend.inc"
#endif

#endif
