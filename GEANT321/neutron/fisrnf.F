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
C*********************************************************************
      FUNCTION FISRNF(A,B)
C*********************************************************************
C Sample secondary fission neutron energy from Watt spectrum
C taken from ORNL/TM-7631
C CZ 3/11/92
      DIMENSION RNDM(3)
C
      CALL GRNDM(RNDM,3)
      Z=SQRT(-ALOG(RNDM(1)))
      S=6.28319*RNDM(2)
      ALOGR3=ALOG(RNDM(3))
      X=SQRT(A*B)/2.
      E1=A*((Z*COS(S)+X)**2-ALOGR3)
C--  E2=A*((Z*SIN(S)+X)**2-ALOGR3)
C distribution of E1 and E2 are identical
      FISRNF = E1
      RETURN
      END
