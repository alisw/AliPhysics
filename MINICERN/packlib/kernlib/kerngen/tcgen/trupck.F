*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/02/15 17:49:56  mclareni
* Kernlib
*
*
#include "kerngen/pilot.h"
      SUBROUTINE TRUPCK (U,S,M)
C
C CERN PROGLIB# F112    TRUPCK          .VERSION KERNFOR  2.08  741218
C ORIG. 18/12/74 WH
C
      COMMON /SLATE/ I,IM,IS,IU,IV,IH,M2,DUM(33)
      DIMENSION U(*),S(*)
C
      M2 = M*M
      IS = M2
      IU = (M2+M)/2
      I = M - 1
C
   10 IM = I*M
   20 S(IS) = U(IU)
      IS = IS - 1
      IU = IU - 1
      IF (IS.GT.IM)           GO TO 20
      IS = IS - M + I
      I = I - 1
      IF (I.GE.0)             GO TO 10
C
      IS = 1
   40 IV = IS
      IH = IS
   50 IV = IV + M
      IH = IH + 1
      IF (IV.GT.M2)           GO TO 60
      S(IH) = S(IV)
        GO  TO                50
   60 IS = IS + M + 1
      IF (IS.LT.M2)           GO TO 40
C
      RETURN
      END
