*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/04/01 15:02:51  mclareni
* Mathlib gen
*
*
#include "gen/pilot.h"
      SUBROUTINE ORDRE2(A,B,N)
C              IT ORDERS N PAIRS OF NUMBERS (X(I),Y(I))  IN ASCENDING
C              ORDER OF X(I)
C
      DIMENSION A(1000), B(1000)
      IF  (N .LE. 1)   RETURN
      NM1 = N-1
      DO  200  IV= 1, NM1
      IQUIT = 0
      NH = N - IV
      DO  100  IH = 1, NH
      IF  (A(IH) - A(IH+1))  100, 100, 50
   50 SAV = A(IH)
      A(IH) = A(IH+1)
      A(IH+1) = SAV
      SAV = B(IH)
      B(IH) = B(IH+1)
      B(IH+1) = SAV
      IQUIT = 1
  100 CONTINUE
      IF  (IQUIT)  250, 250, 200
  200 CONTINUE
  250 CONTINUE
      RETURN
      END
