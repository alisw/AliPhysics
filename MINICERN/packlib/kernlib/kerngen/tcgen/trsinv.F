*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/02/15 17:49:54  mclareni
* Kernlib
*
*
#include "kerngen/pilot.h"
      SUBROUTINE TRSINV (G,GI,N)
C
C CERN PROGLIB# F112    TRSINV          .VERSION KERNFOR  2.08  741218
C ORIG. 18/12/74 WH
C
      DIMENSION G(2),GI(2)
C
      CALL TRCHLU (G,GI,N)
      CALL TRINV (GI,GI,N)
      CALL TRSMUL (GI,GI,N)
C
      RETURN
      END
