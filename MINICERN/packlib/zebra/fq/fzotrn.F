*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/03/06 10:47:14  mclareni
* Zebra
*
*
#include "zebra/pilot.h"
      SUBROUTINE FZOTRN (MBAT,NWBAT)

C-    Transmit one batch of data to the output buffer,
C-    exchange file format, native data format
C-    service routine to FZOUT
C-    careful : must not be called with  NWBAT=0

#include "zebra/zunit.inc"
#include "zebra/zebq.inc"
#include "zebra/mzca.inc"
#include "zebra/fzcx.inc"
C--------------    End CDE                             --------------
      DIMENSION    MBAT(99), NWBAT(9)
#if (defined(CERNLIB_QTRHOLL))&&(!defined(CERNLIB_A6M))
      DIMENSION    NAMESR(2)
      DATA  NAMESR / 4HFZOT, 4HRN   /
#endif
#if (defined(CERNLIB_QTRHOLL))&&(defined(CERNLIB_A6M))
      DATA  NAMESR / 6HFZOTRN /
#endif
#if !defined(CERNLIB_QTRHOLL)
      CHARACTER    NAMESR*8
      PARAMETER   (NAMESR = 'FZOTRN  ')
#endif


*     quick trace without test on capacity
#include "zebra/qtraceq.inc"

      NWDO  = NWBAT(1)
      NDONE = 0
      IF (IDX(2).NE.0)             GO TO 44

   24 NWAV = L4ENDX - L4CURX
      NW   = MIN (NWDO-NDONE, NWAV)
      CALL UCOPY (MBAT(NDONE+1), LQ(L4CURX), NW)
      L4CURX = L4CURX + NW
      IF (L4CURX.EQ.L4ENDX)        GO TO 41
#include "zebra/qtrace99.inc"
      RETURN

C--                Batch continued

   41 NDONE = NDONE + NW

C--                Start new logical record

   44 CALL FZOREC
      IF (NDONE.LT.NWDO)           GO TO 24
      GO TO 999
      END
*      ==================================================
#include "zebra/qcardl.inc"
