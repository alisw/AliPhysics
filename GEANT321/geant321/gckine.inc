*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:20:31  cernlib
* Geant
*
*
#ifndef CERNLIB_GEANT321_GCKINE_INC
#define CERNLIB_GEANT321_GCKINE_INC
*
*
* gckine.inc
*
#if !defined(CERNLIB_INCLUDE)
      COMMON/GCKINE/IKINE,PKINE(10),ITRA,ISTAK,IVERT,IPART,ITRTYP
     +      ,NAPART(5),AMASS,CHARGE,TLIFE,VERT(3),PVERT(4),IPAOLD
C
#if defined(CERNLIB_TYPE)
#include "geant321/gtkine.inc"
#endif
#endif
*
* gckine.inc
*
#if defined(CERNLIB_INCLUDE)
*
#include "geant321/incbeg.inc"
     +   gckine
#include "geant321/incend.inc"
#endif

#endif
