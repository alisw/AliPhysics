*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/02/15 17:49:49  mclareni
* Kernlib
*
*
#include "kerngen/pilot.h"
      SUBROUTINE LORENF (U,PS,PI,PF)
C
C CERN PROGLIB# U102    LORENF          .VERSION KERNFOR  4.04  821124
C ORIG. 20/08/75 L.PAPE
C
#if !defined(CERNLIB_B48M)
      DOUBLE PRECISION PF4, FN
#endif
      DIMENSION      PS(4),PI(4),PF(4)
C
      IF (PS(4).EQ.U) GO TO 17
      PF4  = (PI(4)*PS(4)-PI(3)*PS(3)-PI(2)*PS(2)-PI(1)*PS(1)) / U
      FN   = (PF4+PI(4)) / (PS(4)+U)
      PF(1)= PI(1) - FN*PS(1)
      PF(2)= PI(2) - FN*PS(2)
      PF(3)= PI(3) - FN*PS(3)
      PF(4)= PF4
      GO TO 18
C
   17 PF(1)= PI(1)
      PF(2)= PI(2)
      PF(3)= PI(3)
      PF(4)= PI(4)
C
   18 CONTINUE
C
      RETURN
C
      END
