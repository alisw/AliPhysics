*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:20:33  cernlib
* Geant
*
*
#ifndef CERNLIB_GEANT321_GCPUSH_INC
#define CERNLIB_GEANT321_GCPUSH_INC
*
*
* gcpush.inc
*
#if !defined(CERNLIB_INCLUDE)
      COMMON/GCPUSH/NCVERT,NCKINE,NCJXYZ,NPVERT,NPKINE,NPJXYZ
      INTEGER       NCVERT,NCKINE,NCJXYZ,NPVERT,NPKINE,NPJXYZ
C
#endif
*
* gcpush.inc
*
#if defined(CERNLIB_INCLUDE)
*
#include "geant321/incbeg.inc"
     +   gcpush
#include "geant321/incend.inc"
#endif

#endif
