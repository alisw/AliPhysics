*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/04/01 15:02:26  mclareni
* Mathlib gen
*
*
#include "gen/pilot.h"
      SUBROUTINE SPLAN1(B1,X,X2,X3,X4,ND)
#include "gen/imp64.inc"
      DIMENSION B1(3)

************************************************************************
*   NORBAS, VERSION: 15.03.1993
************************************************************************
*
*   TEST-ROUTINE  SPLAN1
*   FOR COMPUTING  1-DIMENSIONAL  B-SPLINES OF DEGREE  K = 1
*   BY MEANS OF ANALYTIC FORMULAE
*
************************************************************************

      DO 1 I = 1,3
    1 B1(I)=0
      IF(IABS(ND) .GT. 1) RETURN
      IF(X .LT. X2 .OR. X .GT. X4) RETURN
      IF(ND .EQ. 0) THEN
       IF(X .LT. X3) THEN
        B1(1)=(-X+X3)/X3
        B1(2)=X/X3
        B1(3)=0
       ELSE
        B1(1)=0
        B1(2)=(X-1)/(X3-1)
        B1(3)=(-X+X3)/(X3-1)
       END IF
      ELSE IF(ND .EQ. 1) THEN
       IF(X .LT. X3) THEN
        B1(1)=-1/X3
        B1(2)=1/X3
        B1(3)=0
       ELSE
        B1(1)=0
        B1(2)=1/(X3-1)
        B1(3)=-1/(X3-1)
       END IF
      ELSE IF(ND .EQ. -1) THEN
       IF(X .LT. X3) THEN
        B1(1)=-(X/2-X3)*X/X3
        B1(2)=X**2/(2*X3)
        B1(3)=0
       ELSE
        B1(1)=X3/2
        B1(2)=(X**2-2*X+X3)/(2*(X3-1))
        B1(3)=-(X-X3)**2/(2*(X3-1))
       END IF
      END IF
      RETURN
      END



