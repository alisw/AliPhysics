*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:20:35  cernlib
* Geant
*
*
#ifndef CERNLIB_GEANT321_LUDAT3_INC
#define CERNLIB_GEANT321_LUDAT3_INC
*
*
* ludat3.inc
*
#if !defined(CERNLIB_INCLUDE)
      COMMON/LUDAT3/MDCY(500,3),MDME(2000,2),BRAT(2000),KFDP(2000,5)
      INTEGER MDCY,MDME,KFDP
      REAL    BRAT
      SAVE /LUDAT3/
*
#endif
*
* ludat3.inc
*
#if defined(CERNLIB_INCLUDE)
*
#include "geant321/incbeg.inc"
     +    ludat3
#include "geant321/incend.inc"
#endif

#endif
