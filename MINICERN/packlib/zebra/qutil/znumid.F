*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/03/06 10:47:15  mclareni
* Zebra
*
*
#include "zebra/pilot.h"
      SUBROUTINE ZNUMID (ID,NID)

#include "zebra/quest.inc"
C--------------    END CDE                             --------------
      DIMENSION    ID(9), NID(9)


      NAME = ID(1)
      IF (IABS(NAME).GE.10000)     GO TO 22
      CALL VBLANK (IQUEST(26),4)
      CALL USET (NAME,IQUEST(26),1,4)
      CALL UBUNCH (IQUEST(26),NAME,4)
   22 NID(1) = NAME
      RETURN
      END
*      ==================================================
#include "zebra/qcardl.inc"
