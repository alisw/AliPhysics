*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/02/15 17:50:12  mclareni
* Kernlib
*
*
#include "kerngen/pilot.h"
      FUNCTION LVMAX (A,N)
C
C CERN PROGLIB# F121    LVMAX           .VERSION KERNFOR  1.0   720321
C ORIG. 01/07/71
C
      DIMENSION A(*)
C
      IF (N.LE.0)  GO TO 100
      J= 1
      AMAX = A(1)
      DO 9 I= 1,N
      IF (A(I).LE.AMAX)  GO TO 9
      AMAX= A(I)
      J= I
    9 CONTINUE
      LVMAX = J
      RETURN
  100 LVMAX=0
      RETURN
      END
