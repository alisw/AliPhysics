*
* $Id$
*
* $Log$
* Revision 1.2  1996/04/24 17:26:55  mclareni
* Extend the include file cleanup to dzebra, rz and tq, and also add
* dependencies in some cases.
*
* Revision 1.1.1.1  1996/03/06 10:47:24  mclareni
* Zebra
*
*
#include "zebra/pilot.h"
      SUBROUTINE RZINS(IXDIVP,LSUPP,JBIASP,LBANK)
*
************************************************************************
*
*         To read a data structure and make relocation
*
* Called by RZIN
*
*  Author  : J.Zoll/R.Brun/B.Holl
*  Written : 09.05.86
*  Last mod: 20.01.93 RB bug fix in call to MZIOCR
*
************************************************************************
#include "zebra/zmach.inc"
#include "zebra/mqsys.inc"
#include "zebra/mzct.inc"
#include "zebra/mzioc.inc"
#include "zebra/mzcn.inc"
#include "zebra/fzci.inc"
#include "zebra/fzcseg.inc"
#include "zebra/fzcocc.inc"
#include "zebra/rzclun.inc"
      DIMENSION    IXDIVP(9),LSUPP(9),JBIASP(9),IDUM(3)
      EQUIVALENCE (IOPTR,IQUEST(95))
*
*-----------------------------------------------------------------------
*

#if !defined(CERNLIB_FQXISN)
#include "zebra/q_jbyt.inc"
#endif

      IXDIVI = IXDIVP(1)
      JRETCD = 0
      JERROR = 0
      NQOCC  = 0
      NQSEG  = 0
*
*            Read header
*
      CALL RZREAD(NWTABI,3,1,1)
      IF(IQUEST(1).NE.0) GO TO 99
      NWIOI  = 0
      NWUHI  = 0
      NWUHCI = 0
      NWUMXI = 0
      NWTXI  = 0
      NWSEGI = 0
*
*            Get store / divison
*
      CALL MZSDIV(IXDIVI,7)
      IF(JQDIVI.EQ.0) JQDIVI=2
*
      IF(LBANK.NE.0.AND.IOPTR.NE.0)THEN
         NLINK=IQ(KQS+LBANK-3)
         NWK  =NWBKI-10-NLINK
         IF(IQ(KQS+LBANK-1).LT.NWK)GO TO 91
#if !defined(CERNLIB_FQXISN)
         IF(IMODEX.GT.0)THEN
            CALL MZIOCR(LQ(KQS+LBANK-NLINK-1))
            IF(IQUEST(1).LT.0)GO TO 99
            IQUEST(1)=0
         ENDIF
#endif
         CALL RZREAD(IQ(KQS+LBANK+1),NWK,NWTABI+14+NLINK,0)
         GO TO 99
      ENDIF
*
*            Read memory occupation table
*
      LQSYSR(KQT+1)=LSUPP(1)
      CALL FZIMTB
      IF(JRETCD.EQ.3)GO TO 91
      CALL RZREAD(LQ(LQTA+NWTABI),NWTABI,4,1)
      IF(IQUEST(1).NE.0) GO TO 99
*
*            Read the data
*
      LSTA = LQ(LQMTA+3)
      LEND = LQ(LQMTA+4)
#if !defined(CERNLIB_FQXISN)
      IF(IMODEX.LE.0)GO TO 30
*
*            in eXchange mode
*
      LIN = LSTA
      NWR = NWTABI+4
   10 CONTINUE
      CALL RZREAD(LQ(KQS+LIN),1,NWR,1)
      IF(IQUEST(1).NE.0)GO TO 99
      NWR   = NWR+1
      IWD   = LQ(KQS+LIN)
      NST   = JBYT(IWD,1,16)-12
      IF(NST.LT.0)GO TO 20
*
*            True bank
*
      IQLN = LIN
      IQLS = LIN + NST + 1
      IF(IQLS+8.GE.LEND)GO TO 92
      MFO(1) =  1
      MFO(2) =  NST + 2
      MFO(3) =  2
      MFO(4) =  2
      MFO(5) =  5
      MFO(6) =  1
      MFO(7) =  1
      MFO(8) = -1
      JFOEND =  8
      NWFOTT = NST+9
      NWFODN = 0
      CALL RZREAD(LQ(KQS+LIN+1),NST+9,NWR,0)
      IF(IQUEST(1).NE.0)GO TO 99
      NWR    = NWR+NST+9
      IQNIO  = JBYT(IQ(KQS+IQLS),19,4)
      IQNL   = IQ(KQS+IQLS-3)
      IQND   = IQ(KQS+IQLS-1)
      IF(IQNIO+IQNL.NE.NST)GO TO 92
      LIN    = IQLS + IQND + 9
      IF(IQND.GT.0)THEN
         IF(LIN.GT.LEND)GO TO 92
         CALL MZIOCR(LQ(KQS+IQLN))
         IF(IQUEST(1).LT.0)GO TO 99
         IQUEST(1)=0
         CALL RZREAD(IQ(KQS+IQLS+1),IQND,NWR,0)
         IF(IQUEST(1).NE.0)GO TO 99
         NWR = NWR+IQND
      ENDIF
      IF(LIN.LT.LEND)GO TO 10
      GO TO 40
*
*            Short dead region
*
   20 NWD = JBYT(IWD,17,IQDROP-17)
      IF(NWD.EQ.0.OR.NWD.NE.NST+12)GO TO 92
      IF(JBYT(IWD,IQDROP,IQBITW-IQDROP).NE.1)GO TO 92
      LIN = LIN + NWD
      IF(LIN.LT.LEND)GO TO 10
      GO TO 40
#endif
*
*            in native mode
*
   30 NWR = LEND - LSTA
      CALL RZREAD(LQ(KQS+LSTA),NWR,NWTABI+4,0)
      IF(IQUEST(1).NE.0) GO TO 99
 
*
*            Relocation
*
   40 CONTINUE
      CALL FZIREL
      IF(JRETCD.NE.0)GO TO 93
*
*            Link structure and return
*
      JB=JBIASP(1)
      IF(JB.GE.2)THEN
         LSUPP(1)=LENTRI
      ELSE
         LSUPP(1)=LQSYSR(KQT+1)
         CALL ZSHUNT(IXDIVI,LENTRI,LSUPP,JB,1)
      ENDIF
*
      IQUEST(1)  = 0
      IQUEST(11) = IEVFLI
      IQUEST(12) = 0
      IQUEST(13) = LENTRI
      IQUEST(14) = NWBKI
      GO TO 99
*
*            Not enough space
*
 
   91 IQUEST(11)= -2
      IQUEST(1) =  1
      GO TO 99
*
*            Inconsistent bank material
*
   92 IQUEST(11)= -3
      IQUEST(1) =  1
      GO TO 99
*
*            Bad data
*
 
   93 IQUEST(11)= -3
      IQUEST(1) =  1
 
   99 RETURN
      END
