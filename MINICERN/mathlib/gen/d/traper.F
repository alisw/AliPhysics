*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/04/01 15:02:15  mclareni
* Mathlib gen
*
*
#include "gen/pilot.h"
      SUBROUTINE TRAPER (X,Y,E,N,AA,BB,RE,SD)
      DIMENSION X(N),Y(N),E(N)
      RE=0.
      SD=0.
      IF(AA .EQ. BB .OR. N .LT. 2) RETURN
      A=AA
      B=BB
      IF(AA .LT. BB) GO TO 1
      A=BB
      B=AA
    1 DO 40 I1 = 2,N
      I=I1
      IF(X(I) .GE. A) GO TO 41
   40 CONTINUE
   41 DO 42 J1 = 2,N
      J=J1
      IF(X(J) .GE. B) GO TO 43
   42 CONTINUE
   43 WI1=(X(I)-A)**2/(X(I)-X(I-1))
      WI=(1.0+(A-X(I-1))/(X(I)-X(I-1)))*(X(I)-A)
      WJ1=(1.0+(X(J)-B)/(X(J)-X(J-1)))*(B-X(J-1))
      WJ=(B-X(J-1))**2/(X(J)-X(J-1))
      IF(I .NE. J) GO TO 2
      WI1=WI1+WJ1+X(I-1)-X(I)
      WI=WI+WJ+X(I-1)-X(I)
      WJ1=0.
      WJ=0.
      GO TO 10
    2 IF(I .NE. J-1) GO TO 3
      WI=WI+WJ1
      WJ1=0.
      GO TO 10
    3 WI=WI+X(I+1)-X(I)
      WJ1=WJ1+X(J-1)-X(J-2)
      IF(I .EQ. J-2) GO TO 10
      LI=I+1
      LJ=J-2
      DO 4 L = LI,LJ
      RE=RE+(X(L+1)-X(L-1))*Y(L)
    4 SD=SD+((X(L+1)-X(L-1))*E(L))**2
   10 RE=RE+WI1*Y(I-1)+WI*Y(I)+WJ1*Y(J-1)+WJ*Y(J)
      SD=SD+(WI1*E(I-1))**2+(WI*E(I))**2+(WJ1*E(J-1))**2+(WJ*E(J))**2
      RE=0.5*RE
      SD=0.5*SQRT(SD)
      IF(AA .GT. BB) RE=-RE
      RETURN
      END
