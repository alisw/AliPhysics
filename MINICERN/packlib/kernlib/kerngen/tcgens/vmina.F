*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/02/15 17:50:18  mclareni
* Kernlib
*
*
#include "kerngen/pilot.h"
      FUNCTION VMINA (A,N)
C
C CERN PROGLIB# F121    VMINA           .VERSION KERNFOR  1.0   710701
C ORIG. 01/07/71
C
      DIMENSION A(*)
C
      IF (N.LE.0)  GO TO 100
      FMIN= ABS (A(1))
      DO 9 I= 1,N
      FMIN =  MIN (ABS (A(I)), FMIN)
    9 CONTINUE
      VMINA= FMIN
      RETURN
  100 VMINA = 0.
      RETURN
      END
