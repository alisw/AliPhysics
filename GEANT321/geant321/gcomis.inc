*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:20:34  cernlib
* Geant
*
*
#ifndef CERNLIB_GEANT321_GCOMIS_INC
#define CERNLIB_GEANT321_GCOMIS_INC
*
*
* gcomis.inc
*
#if !defined(CERNLIB_INCLUDE)
      COMMON/GCOMIS/JUINIT,JUGEOM,JUKINE,JUSTEP,JUOUT,JULAST
      DIMENSION JPCOMS(6)
      EQUIVALENCE (JPCOMS,JUINIT)
*
#endif
*
* gcomis.inc
*
#if defined(CERNLIB_INCLUDE)
*
#include "geant321/incbeg.inc"
     +   gcomis
#include "geant321/incend.inc"
#endif

#endif
