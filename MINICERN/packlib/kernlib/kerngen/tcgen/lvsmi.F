*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/02/15 17:49:52  mclareni
* Kernlib
*
*
#include "kerngen/pilot.h"
      FUNCTION LVSMI (A,NDO,INC)
C
C CERN PROGLIB# F121    LVSMI           .VERSION KERNFOR  4.21  890323
C ORIG. 01/03/89 FCA
C
C-    Find Minimum in scattered vector, single precision

      DIMENSION A(*)
C
      HMIN = A(1)
      LMIN = 1
      DO 10 J = INC+1, (NDO-1)*INC+1, INC
        IF(HMIN.GT.A(J)) THEN
          HMIN = A(J)
          LMIN = J
        ENDIF
   10 CONTINUE
      LVSMI = LMIN
      END
