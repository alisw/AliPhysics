*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/03/06 10:47:14  mclareni
* Zebra
*
*
#include "zebra/pilot.h"
#if !defined(CERNLIB_FQXISN)
      SUBROUTINE FZOTRX (MBAT,NWBAT)

C-    Transmit one batch of data to the output buffer,
C-    exchange file format, exchange data format
C-    service routine to FZOUT
C-    careful : must not be called with  NWBAT=0

#include "zebra/zunit.inc"
#include "zebra/mqsys.inc"
#include "zebra/eqlqf.inc"
#include "zebra/mzioc.inc"
#include "zebra/fzcx.inc"
C--------------    End CDE                             --------------
      DIMENSION    MBAT(99), NWBAT(9)
#if (defined(CERNLIB_QTRHOLL))&&(!defined(CERNLIB_A6M))
      DIMENSION    NAMESR(2)
      DATA  NAMESR / 4HFZOT, 4HRX   /
#endif
#if (defined(CERNLIB_QTRHOLL))&&(defined(CERNLIB_A6M))
      DATA  NAMESR / 6HFZOTRX /
#endif
#if !defined(CERNLIB_QTRHOLL)
      CHARACTER    NAMESR*8
      PARAMETER   (NAMESR = 'FZOTRX  ')
#endif


*     quick trace without test on capacity
#include "zebra/qtraceq.inc"

      NWFOTT = NWBAT(1)
      NWFORE = NWFOTT
      NWFODN = 0
      IF (IDX(2).NE.0)             GO TO 44

   24 CALL FZOCV (MBAT,LQ(L4CURX))
      L4CURX = L4CURX + IQUEST(1)

   27 IF (NWFOAV.LE.0)             GO TO 41
  991 IF (IFOCON(1).NE.0)          GO TO 71
#include "zebra/qtrace99.inc"
      RETURN

C--                Batch continued

   41 CONTINUE

C--                Start new logical record

   44 CALL FZOREC
      NWFOAV = L4ENDX - L4CURX
      IF (NWFORE.NE.0)             GO TO 24
      GO TO 991

C------            Conversion Problems

   71 IQ(KQSP+LQFX+24) = IQ(KQSP+LQFX+24) + 1
      IF (IFOCON(1).GT.0)  IFOCON(2) = IFOCON(2)
     +                        + (LOCF(MBAT(1)) - LQSTOR)
#if defined(CERNLIB_QPRINT)
      IF (LOGLVX.GE.1)  WRITE (IQLOG,9071) LUNX,IFOCON
 9071 FORMAT (1X/' FZOTRX.  LUN=',I4,', Conversion problem: ',
#endif
#if (defined(CERNLIB_QPRINT))&&(!defined(CERNLIB_HEX))&&(defined(CERNLIB_B60M))
     F 4O23)
#endif
#if (defined(CERNLIB_QPRINT))&&(!defined(CERNLIB_HEX))&&(!defined(CERNLIB_B60M))
     F 4O17)
#endif
#if (defined(CERNLIB_QPRINT))&&(defined(CERNLIB_HEX))
     F 4Z17)
#endif
      GO TO 999
      END
*      ==================================================
#include "zebra/qcardl.inc"
#endif
