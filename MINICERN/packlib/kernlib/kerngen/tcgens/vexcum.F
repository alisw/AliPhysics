*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/02/15 17:50:17  mclareni
* Kernlib
*
*
#include "kerngen/pilot.h"
      SUBROUTINE VEXCUM  (A,EX,N)
C
C CERN PROGLIB# F121    VEXCUM          .VERSION KERNFOR  1.0   710701
C ORIG. 01/07/71
C
      DIMENSION A(*), EX(3)
C
      IF (N.LE.0) RETURN
C
      DO 19 J=1,N
      EX(1)= MIN(EX(1),A(J))
      EX(2)= MAX(EX(2),A(J))
   19 EX(3)= EX(3) + A(J)
      RETURN
      END
