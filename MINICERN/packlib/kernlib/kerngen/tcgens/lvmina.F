*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/02/15 17:50:13  mclareni
* Kernlib
*
*
#include "kerngen/pilot.h"
      FUNCTION LVMINA (A,N)
C
C CERN PROGLIB# F121    LVMINA          .VERSION KERNFOR  1.0   710701
C ORIG. 01/07/71
C
      DIMENSION A(*)
C
      IF (N.LE.0)  GO TO 100
      J= 1
      AMIN= ABS (A(1))
      DO 9 I= 1,N
      CA= ABS (A(I))
      IF (AMIN.LE.CA)  GO TO 9
      AMIN= CA
      J= I
    9 CONTINUE
      LVMINA = J
      RETURN
  100 LVMINA = 0
      RETURN
      END
