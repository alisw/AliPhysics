*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/02/15 17:49:52  mclareni
* Kernlib
*
*
#include "kerngen/pilot.h"
      SUBROUTINE MXUTY(A,I)
C
C CERN PROGLIB# F110    MXUTY           .VERSION KERNFOR  1.0   650809
C ORIG. 01/01/64 RKB
C
      DIMENSION   A(*)
C
      IF(I.EQ.0) GO TO 100
      J=1
      IF(I.EQ.1) GO TO 2
      DO 1 K=2,I
      A(J)=1.
      J=J+1
      DO 1 L=1,I
      A(J)=0.
    1 J=J+1
    2 A(J)=1.
  100 RETURN
      END
