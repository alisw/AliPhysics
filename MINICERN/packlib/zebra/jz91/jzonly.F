*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/03/06 10:47:17  mclareni
* Zebra
*
*
#include "zebra/pilot.h"
      SUBROUTINE JZONLY (IDV)

C-    PURGE REDUNDANT PROCESSORS TITLES

#include "zebra/mqsysh.inc"
#include "zebra/jzuc.inc"
#include "zebra/jzc.inc"
C--------------    END CDE                             -----------------  ------
      DIMENSION    IDV(99)

#include "zebra/qstorjz.inc"

      NID = IDV(1)
      IFL = 0
      L   = LQJZ - 1

   21 L = LQ(KQS+L)
      IF (L.NE.0)  THEN
          IF (IUCOMP(IQ(KQS+L+1),IDV(2),NID).NE.0)   GO TO 21
          CALL MZDROP (IXSTJZ,L, '.')
          GO TO 21
        ENDIF

#if defined(CERNLIB_QDEBUG)
      IF (IFL.EQ.0)  THEN
          IFL = 1
          L   = LQJZ - 2
          GO TO 21
        ENDIF
#endif
      RETURN
      END
*      ==================================================
#include "zebra/qcardl.inc"
