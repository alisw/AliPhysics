*
* $Id$
*
* $Log$
* Revision 1.2  1996/04/18 16:11:27  mclareni
* Incorporate changes from J.Zoll for version 3.77
*
* Revision 1.1.1.1  1996/03/06 10:47:17  mclareni
* Zebra
*
*
#include "zebra/pilot.h"
      SUBROUTINE MZEND

C-    Memory manager termination,  user called

#include "zebra/zbcd.inc"
#include "zebra/zmach.inc"
#include "zebra/zstate.inc"
#include "zebra/zunit.inc"
#include "zebra/mqsys.inc"
C--------------    End CDE                             --------------

#include "zebra/q_shiftl.inc"

      NSTORS = NQSTOR + 1
      IF (NQSTOR.LT.0)       RETURN
      IF (NQSTOR.GE.16)      RETURN

      WRITE (IQLOG,9001) NSTORS
 9001 FORMAT (1X/' MZEND.   Usage statistics for',I3,' dynamic stores.')

      JSTO = 0
   21 JDIV = ISHFTL (JSTO,26)

      CALL MZSDIV (JDIV,-7)
      CALL MZMAPS

      JSTO = JSTO + 1
      IF (JSTO.LE.NQSTOR)          GO TO 21
      RETURN
      END
*      ==================================================
#include "zebra/qcardl.inc"
