*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/03/06 10:47:17  mclareni
* Zebra
*
*
#include "zebra/pilot.h"
      SUBROUTINE JZLIFT

C-    Initialize processor SV bank for support variables

C-    Parameters in /JZC/ :
C-    IAM(1)  IAMID  = processor ID in A4
C-        2   IAFLDW = down-call flag
C-             0  no further down transfer
C-             1  with down transfer, level JQMLEV must not be reached
C-        3   IANAN  = number of processor constants
C-        4   IANCR  = number of conditions to be recorded
C-        5   IANLSV = number of wsp links to be saved
C-        6   IANDSV = number of wsp data words to be saved

#include "zebra/zstate.inc"
#include "zebra/zunit.inc"
#include "zebra/mqsysh.inc"
#include "zebra/jzuc.inc"
#include "zebra/jzc.inc"
C--------------    END CDE                             -----------------  ------
      DIMENSION    IDPR(7)
      EQUIVALENCE (IDPR(1),IAMID)

      DIMENSION    MMSV(5)
#if (defined(CERNLIB_QTRHOLL))&&(!defined(CERNLIB_A6M))
      DIMENSION    NAMESR(2)
      DATA  NAMESR / 4HJZLI, 4HFT   /
#endif
#if (defined(CERNLIB_QTRHOLL))&&(defined(CERNLIB_A6M))
      DATA  NAMESR / 6HJZLIFT /
#endif
#if !defined(CERNLIB_QTRHOLL)
      CHARACTER    NAMESR*8
      PARAMETER   (NAMESR = 'JZLIFT  ')
#endif
      DATA  MMSV   / 4HJZSV, -7, 0, -7, 0 /

#include "zebra/q_jbit.inc"

#include "zebra/qtraceq.inc"

C----              Create SV bank

      IANCR   = MAX (IANCR,1)
      MMSV(4) = JQNACC + IANCR + IANAN + IANDSV + 2
      IANFL   = 0
#if defined(CERNLIB_QDEBUG)
      LTFL   = LZFIND (IXSTJZ,LQ(KQS+LQJZ-2), IAMID,1)
      IF (LTFL.NE.0)  IANFL = IQ(KQS+LTFL-1) - 1
      MMSV(4) = MMSV(4) + IANFL
#endif
      MMSV(2) = IANLSV + 3
      CALL MZLIFT (IXDVJZ,LSV,LQJZ,-3,MMSV,0)
      LCR = LSV + JQNACC
      LAN = LCR + IANCR + 1
#if defined(CERNLIB_QDEBUG)
      LFL = LAN + IANAN + 1 + IANDSV
#endif
      CALL SBYT (JQLEV+1, IQ(KQS+LSV),1,8)
      IQ(KQS+LSV+1) = IAMID
      IQ(KQS+LSV+3) = IANLSV
      IQ(KQS+LSV+4) = IANDSV

      IQ(KQS+LCR) = IANCR
      IQ(KQS+LAN) = IANAN
      NANTL    = IANAN

C----              Find and copy JQAN title, if any

      LTAN = LZFIND (IXSTJZ,LQ(KQS+LQJZ-1), IAMID,1)
      IF (LTAN.NE.0)  THEN
          NANTL = IQ(KQS+LTAN-1) - 1
          N     = MIN (NANTL,IANAN)
          CALL UCOPY (IQ(KQS+LTAN+2),IQ(KQS+LAN+1),N)
          CALL SBIT1 (IQ(KQS+LSV),17)
          CALL MZDROP (IXSTJZ,LTAN, '.')
        ENDIF

C----              Find and copy JQFL flag titles, if any

#if defined(CERNLIB_QDEBUG)
      IF (LTFL.NE.0)  THEN
          LTFL = LZFIND (IXSTJZ,LQ(KQS+LQJZ-2), IAMID,1)
          IF (IANFL.NE.0)
     +      CALL UCOPY (IQ(KQS+LTFL+2),IQ(KQS+LFL+1),IANFL)
          IQ(KQS+LFL) = IANFL
          CALL MZDROP (IXSTJZ,LTFL, '.')
        ENDIF
#endif

C----              Print and check discrepancies

      LEV = JQLLEV
      IF (IANAN.NE.NANTL)  LEV = LEV + 1

      IF (LEV.GE.2)  WRITE (IQLOG,9042) IDPR
 9042 FORMAT (/' JZLIFT.  Init of "',A4,
     F ' with  IFDW,NAN,NCD,NLSV,NDSV,NFL =',6I4)

      IF (IANAN.EQ.NANTL)          GO TO 47

      IF (LEV.GE.2)  WRITE (IQLOG,9043) NANTL
 9043 FORMAT (10X,'!!! NAN from title =',I8,' !!!')

C----              Exit

   47 IQUEST(1) = -7
      IQUEST(2) = LSV
      IF (JBIT(IQ(KQS+LSV),17).EQ.0)  GO TO 999
      IF (JQLLEV.GE.2)  WRITE (IQLOG,9048)
 9048 FORMAT (10X,'with title.')
      IQUEST(1) = 0
#include "zebra/qtrace99.inc"
      RETURN
      END
*      ==================================================
#include "zebra/qcardl.inc"
