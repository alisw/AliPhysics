*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/03/06 10:47:14  mclareni
* Zebra
*
*
#include "zebra/pilot.h"
      SUBROUTINE ZABEND

C-    DEFAULT FOR THE USER ROUTINE TO HANDLE ABNORMAL PROGRAM END

#include "zebra/zstate.inc"
#include "zebra/zunit.inc"

      CALL ZPOSTM ('TCWM.')
      IF (NQPHAS.GT.0)             GO TO 12
      CALL ABEND
   12 NQPHAS = -2
      CALL ZEND
      RETURN
      END
*      ==================================================
#include "zebra/qcardl.inc"
