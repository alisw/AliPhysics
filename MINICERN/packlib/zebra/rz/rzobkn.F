*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/03/06 10:47:25  mclareni
* Zebra
*
*
#include "zebra/pilot.h"
      SUBROUTINE RZOBKN
*
************************************************************************
*
*      RZ internal routine to write data according to internal table
*      previously built by FZOTAB
*
* Called by RZOUT
*
*  Author  : J.Zoll,R.Brun,B.Holl
*  Written : 03.05.86
*  Last mod: 25.06.93 JDS Add IEVENT to RZWRT calls
*
************************************************************************
#include "zebra/mqsys.inc"
#include "zebra/mzct.inc"
#include "zebra/mzcn.inc"
#include "zebra/mzioc.inc"
#include "zebra/fzcx.inc"
#include "zebra/rzclun.inc"
*
*-----------------------------------------------------------------------
*
*     New event
*
      IEVENT = 0
*
*         Write NWTABX,NWBKX,LENTRX
*
#if !defined(CERNLIB_FQXISN)
      IF(IMODEX.GT.0)THEN
         MFO(1) = 1
         MFO(2) = -1
         JFOEND = 2
      ENDIF
#endif
      CALL RZWRT(NWTABX,3,0,IEVENT)
      IF(IQUEST(1).NE.0) GO TO 99
*
*     Event continued
*
      IEVENT = 1
*
*         Write memory description table
*
      IF(NWTABX.GT.0)THEN
         CALL RZWRT(LQ(LQTA),NWTABX,0,IEVENT)
         IF(IQUEST(1).NE.0) GO TO 99
      ENDIF
*
*         Write bank material
*
      LTB = LQTA
  10  L   = LQ(LTB)
      LE  = LQ(LTB+1)
#if !defined(CERNLIB_FQXISN)
      IF(IMODEX.LE.0)GO TO 50
*
*        in eXchange format
*
  20  IQLN= L
      IWD = LQ(KQS+L)
      NST = JBYT(IWD,1,16) - 12
      IF(NST.LT.0)GO TO 40
*
*         True bank
*
      IQLS  = L + NST + 1
      IQNIO = JBYT(IQ(KQS+IQLS),19,4)
      IQNL  = IQ(KQS+IQLS-3)
      IQND  = IQ(KQS+IQLS-1)
      IQNX  = IQLS + IQND + 9
*
*         First word, I/O words, links, next-link, up-link
*
      N = IQNIO + IQNL
      MFO(1) = 1
      MFO(2) = N + 3
      MFO(3) = 2
      MFO(4) = 2
      MFO(5) = 5
      MFO(6) = 1
      MFO(7) = 1
      MFO(8) = -1
      JFOEND = 8
      CALL RZWRT(LQ(KQS+L),N+10,0,IEVENT)
      IF(IQUEST(1).NE.0) GO TO 99
*
*         Data words
*
      IF( IQND .NE. 0 ) THEN
         CALL MZIOCR(LQ(KQS+IQLN))
         IF(IQUEST(1).LT.0)GO TO 99
         IQUEST(1)=0
         CALL RZWRT(IQ(KQS+IQLS+1),IQND,0,IEVENT)
         IF(IQUEST(1).NE.0) GO TO 99
      ENDIF
   30 L = IQNX
      IF(L.LT.LE) GO TO 20
      LTB = LTB + 2
      IF(LTB.LT.LQTE) GO TO 10
      GO TO 99
*
*         Short dead region
*
   40 NWD = NST + 12
      IQNX= L + NWD
      GO TO 30
#endif
*
*         in Native mode
*
   50 N = LE - L
      IF(N.GT.0)THEN
         CALL RZWRT(LQ(KQS+L),N,0,IEVENT)
         IF(IQUEST(1).NE.0) GO TO 99
      ENDIF
      LTB = LTB + 2
      IF(LTB.LT.LQTE)GO TO 10
*
   99 RETURN
      END
