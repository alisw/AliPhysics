*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/03/06 10:47:21  mclareni
* Zebra
*
*
#include "zebra/pilot.h"
      SUBROUTINE MZIOBK (MMBK,NWBK,CHFORM)

C-    Analyse I/O descriptor,  user called
C-    store the resulting I/O characteristic into the user's
C-    bank description vector MMBK of length NWBK, starting at word 5

#include "zebra/zkrakc.inc"
C--------------    END CDE                             --------------
      DIMENSION    MMBK(99), NWBK(9)
      CHARACTER    CHFORM*(*)


      NW = NWBK(1) - 4
      IQCETK(121) = MMBK(1)
      CALL MZIOCH (MMBK(5),NW,CHFORM)
      RETURN
      END
*      ==================================================
#include "zebra/qcardl.inc"
