*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:20:33  cernlib
* Geant
*
*
#ifndef CERNLIB_GEANT321_GCUNIT_INC
#define CERNLIB_GEANT321_GCUNIT_INC
*
*
* gcunit.inc
*
#if !defined(CERNLIB_INCLUDE)
      COMMON/GCUNIT/LIN,LOUT,NUNITS,LUNITS(5)
      INTEGER LIN,LOUT,NUNITS,LUNITS
      COMMON/GCMAIL/CHMAIL
      CHARACTER*132 CHMAIL
C
#endif
*
* gcunit.inc
*
#if defined(CERNLIB_INCLUDE)
*
#include "geant321/incbeg.inc"
     +   gcunit
#include "geant321/incend.inc"
#endif

#endif
