*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:20:33  cernlib
* Geant
*
*
#ifndef CERNLIB_GEANT321_GCSETS_INC
#define CERNLIB_GEANT321_GCSETS_INC
*
*
* gcsets.inc
*
#if !defined(CERNLIB_INCLUDE)
      COMMON/GCSETS/IHSET,IHDET,ISET,IDET,IDTYPE,NVNAME,NUMBV(20)
C
#if defined(CERNLIB_TYPE)
#include "geant321/gtsets.inc"
#endif
#endif
*
* gcsets.inc
*
#if defined(CERNLIB_INCLUDE)
*
#include "geant321/incbeg.inc"
     +   gcsets
#include "geant321/incend.inc"
#endif

#endif
