*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:20:31  cernlib
* Geant
*
*
#ifndef CERNLIB_GEANT321_GCURSB_INC
#define CERNLIB_GEANT321_GCURSB_INC
*
*
* gcursb.inc
*
#if !defined(CERNLIB_INCLUDE)
      COMMON/GCURSB/NUMNDS,IADDI,NUMND2,NNPAR,IISELT
      COMMON/GCURSC/MOMO
      CHARACTER*4 MOMO
*
#if defined(CERNLIB_TYPE)
#include "geant321/gtursb.inc"
#endif
#endif
*
* gcursb.inc
*
#if defined(CERNLIB_INCLUDE)
*
#include "geant321/incbeg.inc"
     +   gcursb
#include "geant321/incend.inc"
#endif

#endif
