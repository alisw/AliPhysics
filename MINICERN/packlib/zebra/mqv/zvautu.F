*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/03/06 10:47:13  mclareni
* Zebra
*
*
#include "zebra/pilot.h"
      SUBROUTINE ZVAUTU

C-    Push the user VID stack

#include "zebra/zvfaut.inc"
C--------------    END CDE                             --------------

      IF (IQVSTA.EQ.0)             RETURN

      IQVREM(1,6) = IQVREM(1,5)
      IQVREM(2,6) = IQVREM(2,5)

      IQVREM(1,5) = IQVREM(1,4)
      IQVREM(2,5) = IQVREM(2,4)

      IQVREM(1,4) = IQVID(1)
      IQVREM(2,4) = IQVID(2)
      RETURN
      END
*      ==================================================
#include "zebra/qcardl.inc"
