*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:21:12  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.40  by  S.Giani
*-- Author :
      SUBROUTINE HIGSEL(ISEL)
C
#include "geant321/s_result.inc"
C
      DIMENSION RNDM(1)
      CALL GRNDM(RNDM,1)
      ISEL=1
      IF(P.LT.25.+RNDM(1)*25.) ISEL=0
      RETURN
      END
