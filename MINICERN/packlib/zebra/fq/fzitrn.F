*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/03/06 10:47:15  mclareni
* Zebra
*
*
#include "zebra/pilot.h"
      SUBROUTINE FZITRN (MBAT,NWBAT)

C-    Transmit one batch of data from the input buffer,
C-    exchange file format, native data format
C-    service routine to FZIN
C-    careful : must not be called with  N4SKII=0 .AND. NWBAT=0

#include "zebra/zunit.inc"
#include "zebra/zebq.inc"
#include "zebra/mzca.inc"
#include "zebra/fzci.inc"
C--------------    End CDE                             --------------
      DIMENSION    MBAT(99), NWBAT(9)
#if (defined(CERNLIB_QTRHOLL))&&(!defined(CERNLIB_A6M))
      DIMENSION    NAMESR(2)
      DATA  NAMESR / 4HFZIT, 4HRN   /
#endif
#if (defined(CERNLIB_QTRHOLL))&&(defined(CERNLIB_A6M))
      DATA  NAMESR / 6HFZITRN /
#endif
#if !defined(CERNLIB_QTRHOLL)
      CHARACTER    NAMESR*8
      PARAMETER   (NAMESR = 'FZITRN  ')
#endif


*     quick trace without test on capacity
#include "zebra/qtraceq.inc"

      NWDO  = NWBAT(1)
      NDONE = 0

   24 NWAV = L4ENDI - L4CURI
      IF (NWAV.EQ.0)               GO TO 41
      IF (N4SKII.NE.0)             GO TO 31

      NW  = MIN (NWDO-NDONE, NWAV)
      CALL UCOPY (LQ(L4CURI), MBAT(NDONE+1), NW)
      NDONE  = NDONE  + NW
      L4CURI = L4CURI + NW
      IF (NDONE.LT.NWDO)           GO TO 41
#include "zebra/qtrace99.inc"
      RETURN

C--                Skip over N4SKII leading data words

   31 NSK    = MIN (N4SKII,NWAV)
      L4CURI = L4CURI + NSK
      N4SKII = N4SKII - NSK
      IF (N4SKII.EQ.0)             GO TO 46

C----              Read next PhR / LR for batch continued

   41 IFLAGI = -1
      CALL FZIREC
      IF (IFLAGI.NE.0)             GO TO 999
   46 IF (NWDO+N4SKII.NE.0)        GO TO 24
      GO TO 999
      END
*      ==================================================
#include "zebra/qcardl.inc"
