*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:20:31  cernlib
* Geant
*
*
#ifndef CERNLIB_GEANT321_GCCUTS_INC
#define CERNLIB_GEANT321_GCCUTS_INC
*
*
* gccuts.inc
*
#if !defined(CERNLIB_INCLUDE)
      COMMON/GCCUTS/CUTGAM,CUTELE,CUTNEU,CUTHAD,CUTMUO,BCUTE,BCUTM
     +             ,DCUTE ,DCUTM ,PPCUTM,TOFMAX,GCUTS(5)
C
#if defined(CERNLIB_TYPE)
#include "geant321/gtcuts.inc"
#endif
#endif
*
* gccuts.inc
*
#if defined(CERNLIB_INCLUDE)
*
#include "geant321/incbeg.inc"
     +   gccuts
#include "geant321/incend.inc"
#endif

#endif
