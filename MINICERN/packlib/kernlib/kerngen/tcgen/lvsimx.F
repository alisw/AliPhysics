*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/02/15 17:49:51  mclareni
* Kernlib
*
*
#include "kerngen/pilot.h"
      FUNCTION LVSIMX (A,NDO,INC)
C
C CERN PROGLIB# F121    LVSIMX          .VERSION KERNFOR  4.21  890323
C ORIG. 01/03/89 FCA
C
C-    Find Maximum in scattered vector, integer

      INTEGER A(*), HMAX
C
      HMAX = A(1)
      LMAX = 1
      DO 10 J = INC+1, (NDO-1)*INC+1, INC
        IF(HMAX.LT.A(J)) THEN
          HMAX = A(J)
          LMAX = J
        ENDIF
   10 CONTINUE
      LVSIMX = LMAX
      END
