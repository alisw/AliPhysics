*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/02/15 17:50:17  mclareni
* Kernlib
*
*
#include "kerngen/pilot.h"
      SUBROUTINE VMATL (H,A,X,K,N)
C
C CERN PROGLIB# F121    VMATL           .VERSION KERNFOR  1.0   710701
C ORIG. 01/07/71
C
      DIMENSION A(*),X(*),H(N,K)
C
      IF (N.LE.0)  RETURN
      DO 9 I= 1,K
      X(I)= 0.
      DO 9 J= 1,N
      X(I)= X(I) + A(J)*H(J,I)
    9 CONTINUE
      RETURN
      END
