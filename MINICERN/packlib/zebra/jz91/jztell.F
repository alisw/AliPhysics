*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/03/06 10:47:17  mclareni
* Zebra
*
*
#include "zebra/pilot.h"
      SUBROUTINE JZTELL (IDV)

C-    Processor condition recording

C     IDV(1)  condition number

#include "zebra/zunit.inc"
#include "zebra/mqsysh.inc"
#include "zebra/jzuc.inc"
#include "zebra/jzc.inc"
C--------------    END CDE                             -----------------  ------
      DIMENSION    IDV(7)

#if defined(CERNLIB_QDEBUG)
#include "zebra/q_jbit.inc"
#endif

#include "zebra/qstorjz.inc"

      ID = IDV(1)

#if defined(CERNLIB_QDEBUG)
      IF (JBIT(JQLLEV,5).NE.0)  WRITE (IQLOG,9024) JQLEV,NQME(1),ID
 9024 FORMAT ( ' =======  JZTELL level',I2,', "',A4,I4)
#endif

      IF (LQSV.EQ.0)               RETURN
      L = LQSV + JQNACC
      J = MIN (ID,IQ(KQS+L))
      J = MAX (J,1)
      IQ(KQS+L+J) = IQ(KQS+L+J) + 1
      RETURN
      END
*      ==================================================
#include "zebra/qcardl.inc"
