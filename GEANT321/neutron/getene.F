*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:21:56  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/04 23/02/95  14.46.01  by  S.Giani
*-- Author :
      SUBROUTINE GETENE(EN,D1,D2,N)
C sample N times secondary energy distribution and
C store in ENE(*)
#include "geant321/mconst.inc"
#include "geant321/mnutrn.inc"
C
      DIMENSION D1(*),D2(*)
C
      DO 10 I=1,N
         CALL SECEGY(EX,D1,EN,D2)
         ENE(I) = EX
  10  CONTINUE
      RETURN
      END
 
