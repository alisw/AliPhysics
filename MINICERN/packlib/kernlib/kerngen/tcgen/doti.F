*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/02/15 17:49:51  mclareni
* Kernlib
*
*
#include "kerngen/pilot.h"
      FUNCTION DOTI (X,Y)
C
C CERN PROGLIB# F116    DOTI            .VERSION KERNFOR  3.03  780429
C ORIG. 01/07/71
C
      DIMENSION X(4),Y(4)
C
      DOTI = X(1)*Y(1) + X(2)*Y(2) + X(3)*Y(3) - X(4)*Y(4)
      RETURN
      END
