*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/03/06 10:47:14  mclareni
* Zebra
*
*
#include "zebra/pilot.h"
      SUBROUTINE ZFATAM (MESSAG)

C-    FATAL PROGRAM TERMINATION WITH MESSAGE

#include "zebra/mzca.inc"
#include "zebra/zstate.inc"
#include "zebra/zunit.inc"
C--------------    END CDE                             --------------
      CHARACTER    MESSAG*(*)
#if defined(CERNLIB_QMVDS)
      SAVE         INIT
#endif
#include "zebra/zfatalch.inc"

#include "zebra/zfatalre.inc"

#if defined(CERNLIB_QPRINT)
      IF (NQERR.GE.3)  GO TO 19
      NQERR = NQERR+1

      WRITE (IQLOG,9001) MESSAG

      IF (IQTYPE.EQ.0)             GO TO 18
      IF (IQTYPE.EQ.IQLOG)         GO TO 18
      WRITE (IQTYPE,9001) MESSAG
 9001 FORMAT (1X/' ZFATAM.  !!!!!  Going to ZFATAL for ',A)

   18 NQERR = NQERR - 1
#endif
   19 CALL ZFATAL

      END
*      ==================================================
#include "zebra/qcardl.inc"
