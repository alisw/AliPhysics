*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/02/15 17:50:35  mclareni
* Kernlib
*
*
#include "kerngen/pilot.h"
      SUBROUTINE XOUTBF (LUN,BV,NB)
C
C CERN PROGLIB# Z203    XOUTBF          .VERSION KERNFOR  3.05  780715
C ORIG. 15/07/78 JZ
C
      DIMENSION    BV(NB)
C
      WRITE (LUN)  BV
      RETURN
      END
