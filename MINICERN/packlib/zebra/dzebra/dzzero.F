*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/03/06 10:47:07  mclareni
* Zebra
*
*
*-----------------------------------------------------------
#include "zebra/pilot.h"
#if (defined(CERNLIB_DEBUGON))&&(defined(CERNLIB_VFORT))
#include "zebra/debugvf1.inc"
#endif
      SUBROUTINE DZZERO (IXSTOR,LBANK)
#include "zebra/mqsys.inc"
#include "zebra/mzcn.inc"
#include "zebra/mzioc.inc"

      CHARACTER        CLZERO*4
      DOUBLE PRECISION DLZERO
      INTEGER          HLZERO , ILZERO
      REAL             RLZERO
      SAVE CLZERO,DLZERO,HLZERO,ILZERO,RLZERO,IFIRST
      DATA CLZERO,DLZERO,ILZERO,RLZERO /'    ',0.0D0,0,0.0/ ,IFIRST /1/

#if (defined(CERNLIB_DEBUGON))&&(defined(CERNLIB_VFORT))
#include "zebra/debugvf2.inc"
#endif




      IF (IFIRST.EQ.1) THEN
          IFIRST=0
          CALL UCTOH(CLZERO,HLZERO,4,4)
      ENDIF


      CALL MZSDIV(IXSTOR,-1)
      CALL MZCHLS(-7,LBANK)
      CALL MZIOCR(LQ(KQS+IQLN))
      IF(IQFOUL.NE.0)                                      GO TO 999

      JDATA  = 0
      JFOCUR = 0


   10 ITYPE  = MFO(JFOCUR+1)
      IF (ITYPE.EQ.7)                                      GO TO 40
      NWSEC = MFO(JFOCUR+2)
      IF   (NWSEC)           20, 30, 60


   20 NWSEC = IQND - JDATA
                                                           GO TO 60


   30 JDATA = JDATA + 1
      IWORD = IQ(KQS+LBANK+JDATA)
      NWSEC = IWORD
                                                           GO TO 50


   40 JDATA = JDATA + 1
      IWORD = IQ(KQS+LBANK+JDATA)
      ITYPE = MOD (IWORD,16)
      NWSEC = IWORD/16

   50 IF (ITYPE.GT.5.OR.ITYPE.LT.0) THEN
          ITYPE      = 0
          NWSEC      = IQND - JDATA
      ELSEIF (NWSEC.LE.0) THEN
          NWSEC      = IQND - JDATA
      ENDIF


   60 IDBLE  = 0
      DO 100 I=JDATA+1,JDATA+NWSEC
          IF (IDBLE.EQ.1)           THEN
              IDBLE  = 0
          ELSEIF (ITYPE.EQ.3)  THEN
               Q(KQS+LBANK+I) = RLZERO
          ELSEIF (ITYPE.EQ.4)  THEN
              CALL UCOPY(DLZERO,Q(KQS+LBANK+I),2)
              IDBLE  = 1
          ELSEIF (ITYPE.EQ.5)  THEN
              IQ(KQS+LBANK+I) = HLZERO
          ELSE
              IQ(KQS+LBANK+I) = ILZERO
          ENDIF

  100 CONTINUE


      JDATA = JDATA + NWSEC

      IF (JDATA.GE.IQND)                                   GO TO 999

      IF (JDATA.LT.IQND)   THEN
          JFOCUR = JFOCUR + 2
          IF (JFOCUR.LT.JFOEND)                            GO TO 10
          JFOCUR = JFOREP
                                                           GO TO 10
      ENDIF


  999 RETURN
      END
