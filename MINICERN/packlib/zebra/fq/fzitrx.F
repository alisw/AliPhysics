*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/03/06 10:47:15  mclareni
* Zebra
*
*
#include "zebra/pilot.h"
#if !defined(CERNLIB_FQXISN)
      SUBROUTINE FZITRX (MBAT,NWBAT)

C-    Transmit one batch of data from the input buffer,
C-    exchange file format, exchange data format
C-    service routine to FZIN
C-    careful : must not be called with  N4SKII=0 .AND. NWBAT=0

#include "zebra/zunit.inc"
#include "zebra/mqsys.inc"
#include "zebra/eqlqf.inc"
#include "zebra/mzioc.inc"
#include "zebra/fzci.inc"
C--------------    End CDE                             --------------
      DIMENSION    MBAT(99), NWBAT(9)
#if (defined(CERNLIB_QTRHOLL))&&(!defined(CERNLIB_A6M))
      DIMENSION    NAMESR(2)
      DATA  NAMESR / 4HFZIT, 4HRX   /
#endif
#if (defined(CERNLIB_QTRHOLL))&&(defined(CERNLIB_A6M))
      DATA  NAMESR / 6HFZITRX /
#endif
#if !defined(CERNLIB_QTRHOLL)
      CHARACTER    NAMESR*8
      PARAMETER   (NAMESR = 'FZITRX  ')
#endif


*     quick trace without test on capacity
#include "zebra/qtraceq.inc"

      NWFOTT = NWBAT(1)
      NWFODN = 0

   24 NWFOAV = L4ENDI - L4CURI
      IF (NWFOAV.EQ.0)             GO TO 41
      IF (N4SKII.NE.0)             GO TO 31

      CALL FZICV (LQ(L4CURI),MBAT)
      L4CURI = L4CURI + IQUEST(1)

      IF (IFOCON(1).NE.0)          GO TO 71
   27 IF (NWFORE.NE.0)             GO TO 41
#include "zebra/qtrace99.inc"
      RETURN

C--                Skip over N4SKII leading data words

   31 NSK    = MIN (N4SKII,NWFOAV)
      L4CURI = L4CURI + NSK
      N4SKII = N4SKII - NSK
      IF (N4SKII.EQ.0)             GO TO 46

C----              Read next PhR / LR for batch continued

   41 IFLAGI = -1
      CALL FZIREC
      IF (IFLAGI.NE.0)             GO TO 999
   46 IF (NWFOTT+N4SKII.NE.0)      GO TO 24
      GO TO 999

C------            Conversion Problems

   71 IQ(KQSP+LQFI+24) = IQ(KQSP+LQFI+24) + 1
#if defined(CERNLIB_QPRINT)
      IF (LOGLVI.GE.1)  WRITE (IQLOG,9071) LUNI,IFOCON
 9071 FORMAT (1X/' FZITRX.  LUN=',I4,', Conversion problem: ',
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
      IFOCON(1) = 0
      GO TO 27
      END
*      ==================================================
#include "zebra/qcardl.inc"
#endif
