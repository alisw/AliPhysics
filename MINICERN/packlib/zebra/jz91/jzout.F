*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/03/06 10:47:16  mclareni
* Zebra
*
*
#include "zebra/pilot.h"
      SUBROUTINE JZOUT (CHPA1)

C-    Processor up transfer

C     CHPA1   processor ID in A4

#include "zebra/zstate.inc"
#include "zebra/zunit.inc"
#include "zebra/mqsys.inc"
#include "zebra/jzuc.inc"
#include "zebra/jzc.inc"
C--------------    END CDE                             -----------------  ------
      CHARACTER    CHPA1*4
#if defined(CERNLIB_A4)
      CHARACTER    CHIAM*4
#endif
#if defined(CERNLIB_A8)
      CHARACTER    CHIAM*8
#endif
#if defined(CERNLIB_EQUHOLCH)
      EQUIVALENCE (CHIAM, IAMID)
#endif
#if (defined(CERNLIB_QTRHOLL))&&(!defined(CERNLIB_A6M))
      DIMENSION    NAMESR(2)
      DATA  NAMESR / 4HJZOU, 4HT    /
#endif
#if (defined(CERNLIB_QTRHOLL))&&(defined(CERNLIB_A6M))
      DATA  NAMESR / 6HJZOUT  /
#endif
#if !defined(CERNLIB_QTRHOLL)
      CHARACTER    NAMESR*8
      PARAMETER   (NAMESR = 'JZOUT   ')
#endif

#include "zebra/q_jbit.inc"

#include "zebra/qtraceq.inc"
#include "zebra/qstorjz.inc"

C----              Check matching ID's

      CHIAM  = CHPA1
      IF (CHIAM.NE.'??? ')  THEN
#if !defined(CERNLIB_EQUHOLCH)
          CALL UCTOH (CHIAM, IAMID,4,4)
#endif
          IF (NQME(1).NE.IAMID)        GO TO 91
        ENDIF

C--                Account time

#if defined(CERNLIB_JZTIME)
#include "zebra/jztimout.inc"
#endif

C----              Step level back

      IF (JQLEV.EQ.0)              GO TO 92
      LEVDW = JQLEV
      JQLEV = JQLEV - 1

      J    = LQJZ - JQLEV - 6
      LQDW = LQ(KQS+J)
      LQUP = LQ(KQS+J+1)
      J    = J - JQMLEV
      LQSV = LQ(KQS+J)
      LCD  = LQSV + JQNACC
      LQAN = LCD + IQ(KQS+LCD) + 1

C--                Copy flags

#if defined(CERNLIB_QDEBUG)
      CALL VZERO (JQFLAG,JQMFLW)
      L = LQAN + IQ(KQS+LQAN) + 1 + IQ(KQS+LQSV+4)
      N = IQ(KQS+L)
      IF (N.NE.0)  CALL UCOPY (IQ(KQS+L+1),JQFLAG,N)

      IF (JBIT(JQLLEV,13).NE.0)
     +         WRITE (IQLOG,9024) LEVDW,NQME(1),IQ(KQS+LQSV+1)
 9024 FORMAT ( ' =======  JZOUT  level',I2,', "',A4,'   up to "',A4)
#endif

C---               Work space for upper proc

      NQME(1) = IQ(KQS+LQSV+1)
      IF (JBIT(IQ(KQS+LQDW),15).NE.0)  THEN
          CALL SBIT0 (IQ(KQS+LQDW),15)
          GO TO 999
        ENDIF
      NLINK = IQ(KQS+LQJZ+2*JQLEV+2)
      NDATA = IQ(KQS+LQJZ+2*JQLEV+3)
      CALL MZWORK (IXSTJZ,LQ(KQS+NLINK+1),LQ(KQS+NDATA-1),0)

      N = IQ(KQS+LQSV+3)
      IF (N.NE.0) CALL UCOPY (LQ(KQS+LQSV-N-3),LQ(KQS+NQREF+1),N)
      N = IQ(KQS+LQSV+4)
      IF (N.EQ.0)                  GO TO 999
      L = LQAN + IQ(KQS+LQAN)
      CALL UCOPY (IQ(KQS+L+1),LQ(KQS+NQLINK+1),N)

#include "zebra/qtrace99.inc"
      RETURN

C----              ID mismatch

   92 NQCASE = 1
   91 NQCASE = NQCASE + 1
      NQFATA = 3
      IQUEST(11) = NQME(1)
      IQUEST(12) = IAMID
      IQUEST(13) = JQLEV
#include "zebra/qtofatal.inc"
      END
*      ==================================================
#include "zebra/qcardl.inc"
