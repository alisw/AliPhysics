*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/02/15 17:49:56  mclareni
* Kernlib
*
*
#include "kerngen/pilot.h"
      SUBROUTINE TRPCK (S,U,N)
C
C CERN PROGLIB# F112    TRPCK           .VERSION KERNFOR  2.08  741218
C ORIG. 18/12/74 WH
C
      COMMON /SLATE/ I,IA,IND,IPIV,IDUM(36)
      DIMENSION S(*),U(*)
C
      IA = 0
      IND = 0
      IPIV = 0
C
      DO 20 I=1,N
      IPIV = IPIV + I
C
   10 IA = IA + 1
      IND = IND + 1
      U(IND) = S(IA)
      IF (IND.LT.IPIV)        GO TO 10
C
   20 IA = IA + N - I
C
      RETURN
      END
