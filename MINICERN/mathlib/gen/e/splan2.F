*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/04/01 15:02:26  mclareni
* Mathlib gen
*
*
#include "gen/pilot.h"
      SUBROUTINE SPLAN2(B2,X,X3,X4,X5,X6,ND)
#include "gen/imp64.inc"
      DIMENSION B2(5)

************************************************************************
*   NORBAS, VERSION: 15.03.1993
************************************************************************
*
*   TEST-ROUTINE  SPLAN2
*   FOR COMPUTING  1-DIMENSIONAL  B-SPLINES OF DEGREE  K = 2
*   BY MEANS OF ANALYTIC FORMULAE
*
************************************************************************

      DO 1 I = 1,5
    1 B2(I)=0
      IF(ND .LT. -1 .OR. ND .GT. 2) RETURN
      IF(X .LT. X3 .OR. X .GT. X6) RETURN
      IF(ND .EQ. 0) THEN
       IF(X3 .LT. X4 .AND. X4 .LT. X5 .AND. X5 .LT. X6) THEN
        IF(X .LT. X4) THEN
         B2(1)=((X-X4)/X4)**2
         B2(2)=X*(-X*(X4+X5)+2*X4*X5)/(X4**2*X5)
         B2(3)=X**2/(X4*X5)
         B2(4)=0
         B2(5)=0
        ELSE IF(X4 .LE. X .AND. X .LT. X5) THEN
         B2(1)=0
         B2(2)=-(X-X5)**2/(X5*(X4-X5))
         B2(3)=(X**2*(X4-X5-1)+2*X*X5-X4*X5)/(X5*(X4-1)*(X4-X5))
         B2(4)=(X-X4)**2/((X4-1)*(X4-X5))
         B2(5)=0
        ELSE
         B2(1)=0
         B2(2)=0
         B2(3)=(X-1)**2/((X4-1)*(X5-1))
         B2(4)=(X-1)*(X*(2-X4-X5)+2*X4*X5-X4-X5)/((X4-1)*(X5-1)**2)
         B2(5)=((X-X5)/(X5-1))**2
        END IF
       ELSE IF(X3 .LT. X4 .AND. X4 .EQ. X5 .AND. X5 .LT. X6) THEN
        IF(X .LT. X4) THEN
         B2(1)=((X-X4)/X4)**2
         B2(2)=-2*(X-X4)*X/X4**2
         B2(3)=(X/X4)**2
         B2(4)=0
         B2(5)=0
        ELSE
         B2(1)=0
         B2(2)=0
         B2(3)=((X-1)/(X4-1))**2
         B2(4)=-2*(X-X4)*(X-1)/(X4-1)**2
         B2(5)=((X-X4)/(X4-1))**2
        END IF
       END IF
      ELSE IF(ND .EQ. 1) THEN
       IF(X3 .LT. X4 .AND. X4 .LT. X5 .AND. X5 .LT. X6) THEN
        IF(X .LT. X4) THEN
         B2(1)=2*(X-X4)/X4**2
         B2(2)=-2*(X*(X4+X5)-X4*X5)/(X4**2*X5)
         B2(3)=2*X/(X4*X5)
         B2(4)=0
         B2(5)=0
        ELSE IF(X4 .LE. X .AND. X .LT. X5) THEN
         B2(1)=0
         B2(2)=-2*(X-X5)/((X4-X5)*X5)
         B2(3)=2*(X*(X4-X5-1)+X5)/((X4-X5)*(X4-1)*X5)
         B2(4)=2*(X-X4)/((X4-X5)*(X4-1))
         B2(5)=0
        ELSE
         B2(1)=0
         B2(2)=0
         B2(3)=2*(X-1)/((X4-1)*(X5-1))
         B2(4)=-2*(X*(X4+X5-2)-X4*X5+1)/((X4-1)*(X5-1)**2)
         B2(5)=2*(X-X5)/(X5-1)**2
        END IF
       ELSE IF(X3 .LT. X4 .AND. X4 .EQ. X5 .AND. X5 .LT. X6) THEN
        IF(X .LT. X4) THEN
         B2(1)=2*(X-X4)/X4**2
         B2(2)=-2*(2*X-X4)/X4**2
         B2(3)=2*X/X4**2
         B2(4)=0
         B2(5)=0
        ELSE
         B2(1)=0
         B2(2)=0
         B2(3)=2*(X-1)/(X4-1)**2
         B2(4)=-2*((X-X4)+(X-1))/(X4-1)**2
         B2(5)=2*(X-X4)/(X4-1)**2
        END IF
       END IF
      ELSE IF(ND .EQ. 2) THEN
       IF(X3 .LT. X4 .AND. X4 .LT. X5 .AND. X5 .LT. X6) THEN
        IF(X .LT. X4) THEN
         B2(1)=2/X4**2
         B2(2)=-2*(X4+X5)/(X4**2*X5)
         B2(3)=2/(X4*X5)
         B2(4)=0
         B2(5)=0
        ELSE IF(X4 .LE. X .AND. X .LT. X5) THEN
         B2(1)=0
         B2(2)=-2/((X4-X5)*X5)
         B2(3)=2*(X4-X5-1)/((X4-X5)*(X4-1)*X5)
         B2(4)=2/((X4-X5)*(X4-1))
         B2(5)=0
        ELSE
         B2(1)=0
         B2(2)=0
         B2(3)=2/((X4-1)*(X5-1))
         B2(4)=-2*(X4+X5-2)/((X4-1)*(X5-1)**2)
         B2(5)=2/(X5-1)**2
        END IF
       ELSE IF(X3 .LT. X4 .AND. X4 .EQ. X5 .AND. X5 .LT. X6) THEN
        IF(X .LT. X4) THEN
         B2(1)=2/X4**2
         B2(2)=-4/X4**2
         B2(3)=2/X4**2
         B2(4)=0
         B2(5)=0
        ELSE
         B2(1)=0
         B2(2)=0
         B2(3)=2/(X4-1)**2
         B2(4)=-4/(X4-1)**2
         B2(5)=2/(X4-1)**2
        END IF
       END IF
      ELSE IF(ND .EQ. -1) THEN
       IF(X3 .LT. X4 .AND. X4 .LT. X5 .AND. X5 .LT. X6) THEN
        IF(X .LT. X4) THEN
         B2(1)=(X**2/3-X*X4+X4**2)*X/X4**2
         B2(2)=-((X4+X5)*X-3*X4*X5)*X**2/(3*X4**2*X5)
         B2(3)=X**3/(3*X4*X5)
         B2(4)=0
         B2(5)=0
        ELSE IF(X4 .LE. X .AND. X .LT. X5) THEN
         B2(1)=X4/3
         B2(2)=-(X**3-3*X**2*X5+3*X*X5**2-X4*X5**2)/(3*(X4-X5)*X5)
         B2(3)=(X**3*(X4-X5-1)+3*X**2*X5-3*X*X4*X5+X4**2*X5)/
     1         (3*(X4-X5)*(X4-1)*X5)
         B2(4)=(X-X4)**3/(3*(X4-X5)*(X4-1))
         B2(5)=0
        ELSE
         B2(1)=X4/3
         B2(2)=X5/3
         B2(3)=(X**3-3*X**2+3*X+X4*X5-X4-X5)/(3*(X4-1)*(X5-1))
         B2(4)=(X**3*(-X4-X5+2)+3*X**2*(X4*X5-1)+3*X*(-2*X4*
     1         X5+X4+X5)-X4**2*X5**2+2*X4**2*X5-X4**2+2*X4*X5**2-
     2         X4*X5-X5**2)/(3*(X4-1)*(X5-1)**2)
         B2(5)=(X-X5)**3/(3*(X5-1)**2)
        END IF
       ELSE IF(X3 .LT. X4 .AND. X4 .EQ. X5 .AND. X5 .LT. X6) THEN
        IF(X .LT. X4) THEN
         B2(1)=(X**2/3-X*X4+X4**2)*X/X4**2
         B2(2)=-(2*X-3*X4)*X**2/(3*X4**2)
         B2(3)=X**3/(3*X4**2)
         B2(4)=0
         B2(5)=0
        ELSE
         B2(1)=X4/3
         B2(2)=X4/3
         B2(3)=(X**3-3*X**2+3*X+X4**2-2*X4)/(3*(X4-1)**2)
         B2(4)=-(2*X+X4-3)*(X-X4)**2/(3*(X4-1)**2)
         B2(5)=(X-X4)**3/(3*(X4-1)**2)
        END IF
       END IF
      END IF
      RETURN
      END



