*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/03/06 10:47:15  mclareni
* Zebra
*
*
#include "zebra/pilot.h"
#if !defined(CERNLIB_NODUMMY)
      SUBROUTINE QNEXT

C-    DUMMY VERSION OF THE USER ROUTINE QNEXT

#include "zebra/zunit.inc"
C--------------    END CDE                             --------------


      LUN = IQTYPE
      IF (LUN.NE.0)                GO TO 12
   11 LUN = IQLOG
   12 WRITE (LUN,9001)
      IF (LUN.NE.IQLOG)            GO TO 11
      STOP

 9001 FORMAT ('0!!!!!    DUMMY ROUTINE QNEXT REACHED.')
      END
*      ==================================================
#include "zebra/qcardl.inc"
#endif
