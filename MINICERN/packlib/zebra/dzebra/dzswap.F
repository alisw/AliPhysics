*
* $Id$
*
* $Log$
* Revision 1.2  1996/04/24 17:26:18  mclareni
* Extend the include file cleanup to dzebra, rz and tq, and also add
* dependencies in some cases.
*
* Revision 1.1.1.1  1996/03/06 10:47:07  mclareni
* Zebra
*
*
*----------------------------------------------------------
#include "zebra/pilot.h"
      SUBROUTINE DZSWAP (IXSTOR,LSTR1,LSTR2,CHOPT)


#include "zebra/mqsys.inc"
#include "zebra/mzcn.inc"
#include "zebra/qequ.inc"
#include "zebra/zstate.inc"
#include "zebra/zunit.inc"
#include "zebra/storparq.inc"
#include "zebra/dzc1.inc"

      CHARACTER    *(*) CHOPT
      SAVE NAMESR
#if (defined(CERNLIB_QTRHOLL))&&(!defined(CERNLIB_A6M))
      DIMENSION    NAMESR(2)
      DATA  NAMESR / 4HZSWA, 4HP    /
#endif
#if (defined(CERNLIB_QTRHOLL))&&(defined(CERNLIB_A6M))
      DATA  NAMESR / 6HZSWAP  /
#endif
#if !defined(CERNLIB_QTRHOLL)
      CHARACTER    NAMESR*8
      DATA  NAMESR /  'ZSWAP   ' /
#endif
#include "zebra/q_jbit.inc"
#include "zebra/q_jbyt.inc"

#include "zebra/qtrace.inc"

      IF (IXSTOR.NE.NCHEKQ) CALL MZSDIV (IXSTOR,-1)

      L1    = LSTR1
      L2    = LSTR2

      IFLAG = INDEX(CHOPT,'R')


      IF (LQSTA(KQT+JQDVSY).EQ.LQEND(KQT+JQDVSY)) IFLAG = 0

#if defined(CERNLIB_QDEVZE)
      IF (NQDEVZ.NE.0) WRITE(IQPRNT,
     +'(''0DEVZE ZSWAP,  ENTRY:IXSTOR,LSTR1,LSTR2,IFLAG= '',Z8,3I8)')
     + IXSTOR,L1,L2,IFLAG
#endif

      IF (L1.EQ.0 .OR. L2.EQ.0)                            GO TO 999

      CALL MZCHLS(NCHEKQ,L1)
      IF (IQFOUL.NE.0)                 THEN
          WRITE(IQPRNT,'(''0ZSWAP - BANK1 INVALID'')')
                                                           GO TO 998
      ENDIF

      CALL MZCHLS(NCHEKQ,L2)
      IF (IQFOUL.NE.0)                 THEN
          WRITE(IQPRNT,'(''0ZSWAP - BANK2 INVALID'')')
                                                           GO TO 998
      ENDIF

      LNEXT1 = LQ(KQS+L1)
      IF (LNEXT1.NE.0) CALL MZCHLS(NCHEKQ,LNEXT1)
      IF (IQFOUL.NE.0)                 THEN
          WRITE(IQPRNT,'(''0ZSWAP - Next of LSTR1 invalid'',2I8)')
     +     L1,LNEXT1
                                                           GO TO 998
      ENDIF
      LUP1   = LQLUP(KQS+L1)
      IF (LUP1.NE.0) CALL MZCHLS(NCHEKQ,LUP1)
      IF (IQFOUL.NE.0)                 THEN
          WRITE(IQPRNT,'(''0ZSWAP - Origin of LSTR1 invalid'',2I8)')
     +     L1,LUP1
                                                           GO TO 998
      ENDIF
      LSUP1  = LQLORG(KQS+L1)
      IF (LQ(KQS+LSUP1).NE.L1)         THEN
          WRITE(IQPRNT,'(''0ZSWAP - R link invalid (@r # LSTR1)'',
     +     2I8)') L1,LQ(KQS+LSUP1)
                                                           GO TO 998
      ENDIF

      LNEXT2 = LQ(KQS+L2)
      IF (LNEXT2.NE.0) CALL MZCHLS(NCHEKQ,LNEXT2)
      IF (IQFOUL.NE.0)                 THEN
          WRITE(IQPRNT,'(''0ZSWAP - Next of LSTR2 invalid'',2I8)')
     +     L2,LNEXT2
                                                           GO TO 998
      ENDIF
      LUP2   = LQLUP(KQS+L2)
      IF (LUP2.NE.0) CALL MZCHLS(NCHEKQ,LUP2)
      IF (IQFOUL.NE.0)                 THEN
          WRITE(IQPRNT,'(''0ZSWAP - Origin of LSTR2 invalid'',2I8)')
     +     L2,LUP2
                                                           GO TO 998
      ENDIF
      LSUP2  = LQLORG(KQS+L2)
      IF (LQ(KQS+LSUP2).NE.L2)         THEN
          WRITE(IQPRNT,'(''0ZSWAP - R link invalid (@r # LSTR2)'',
     +     2I8)') L2,LQ(KQS+LSUP2)
                                                           GO TO 998
      ENDIF

#if defined(CERNLIB_QDEVZE)
      IF (NQDEVZ.NE.0) WRITE(IQPRNT,
     +'(''0DEVZE ZSWAP,  Before SWAP ,L1,N,O,P,@P '',5I8,/,
     +  ''       ZSWAP,  Before SWAP ,L2,N,O,P,@P '',5I8)')
     + L1,LQ(KQS+L1),LQLUP(KQS+L1),LQLORG(KQS+L1),
     + LQ(KQS+LQLORG(KQS+L1)),
     + L2,LQ(KQS+L2),LQLUP(KQS+L2),LQLORG(KQS+L2),
     + LQ(KQS+LQLORG(KQS+L2))
#endif

      IF (LNEXT1.EQ.L2) THEN
          LQ(KQS+L1)      = LNEXT2
          LQ(KQS+L2)      = L1
          LQLORG(KQS+L1)  = L2
          LQLORG(KQS+L2)  = LSUP1
          LQ(KQS+LSUP1)   = L2
          IF (LNEXT2.NE.0) LQLORG(KQS+LNEXT2) = L1
      ELSEIF(L1.EQ.LNEXT2) THEN
          LQ(KQS+L1)      = L2
          LQ(KQS+L2)      = LNEXT1
          LQLORG(KQS+L1)  = LSUP2
          LQLORG(KQS+L2)  = L1
          LQ(KQS+LSUP2)   = L1
          IF (LNEXT1.NE.0) LQLORG(KQS+LNEXT1) = L2
      ELSE
          LQ(KQS+L1)      = LNEXT2
          LQ(KQS+L2)      = LNEXT1
          LQLUP(KQS+L1)   = LUP2
          LQLUP(KQS+L2)   = LUP1
          LQLORG(KQS+L1)  = LSUP2
          LQLORG(KQS+L2)  = LSUP1
          LQ(KQS+LSUP1)   = L2
          LQ(KQS+LSUP2)   = L1
          IF (LNEXT1.NE.0) LQLORG(KQS+LNEXT1) = L2
          IF (LNEXT2.NE.0) LQLORG(KQS+LNEXT2) = L1
      ENDIF

      LSTR1  = L2
      LSTR2  = L1

#if defined(CERNLIB_QDEVZE)
      IF (NQDEVZ.NE.0) WRITE(IQPRNT,
     +'(''0DEVZE ZSWAP,  After SWAP ,LSTR1,N,O,P,@P '',5I8,/,
     +  ''       ZSWAP,  After SWAP ,LSTR2,N,O,P,@P '',5I8)')
     + LSTR1,LQ(KQS+LSTR1),LQLUP(KQS+LSTR1),LQLORG(KQS+LSTR1),
     + LQ(KQS+LQLORG(KQS+LSTR1)),
     + LSTR2,LQ(KQS+LSTR2),LQLUP(KQS+LSTR2),LQLORG(KQS+LSTR2),
     + LQ(KQS+LQLORG(KQS+LSTR2))
#endif


      IF (IFLAG.NE.0)      THEN
          LSYSB  = LQSYSS(KQT+MSYLAQ)
          IF(LSYSB.GT.0)                   THEN
              CALL MZCHLS(NCHEKQ,LSYSB)
              IF (IQFOUL.NE.0)                 THEN
                  WRITE(IQPRNT,
     X            '(''0ZSWAP - LINK AREA SYSTEM BANK INVALID'')')
                                                           GO TO 998
              ENDIF


              NWTAB  = IQ(KQS+LSYSB+MLAUSQ)


              LENTRY = LSYSB + KQS + MLAUSQ

              DO 300 IENTRY = 1,(NWTAB-1)/NLAENQ

                  LLAAR1 = IQ(LENTRY+MLAADQ) + KQS
                  LLAARL = IQ(LENTRY+MLALTQ) + KQS
                  NLANS  = IQ(LENTRY+MLANSQ)
                  JTEMP  = JBIT(NLANS,JLATMQ)
                  NTEMP  = NLATMQ*JTEMP
                  NLANS  = JBYT(IQ(LENTRY+MLANSQ),JLANSQ,NLANSQ) - NTEMP
                  IF(JTEMP.EQ.1.AND.LQ(LLAAR1+MLACTQ-1).EQ.0) GO TO 300
                  DO 200 I=LLAAR1+NTEMP+NLANS , LLAARL-1
                      IF (LQ(I).EQ.L1) THEN
                          LQ(I) = L2
                      ELSEIF(LQ(I).EQ.L2) THEN
                          LQ(I) = L1
                      ENDIF
  200             CONTINUE

                  LENTRY = LENTRY + NLAENQ
  300         CONTINUE

          ENDIF



          IDIV1  = 1
          IF (L1.LT.LQEND(KQT+JQDVLL))                     GO TO 410
          IF (L1.GE.LQEND(KQT+NDVMXQ))                     GO TO 998
          IDIV1  = JQDVSY

  410     IF (L1.LT.LQEND(KQT+IDIV1 ))                     GO TO 420
          IDIV1  = IDIV1  + 1
                                                           GO TO 410

  420     IF (L1.LT.LQSTA(KQT+IDIV1 ))                     GO TO 998


          IDIV2  = 1
          IF (L2.LT.LQEND(KQT+JQDVLL))                     GO TO 510
          IF (L2.GE.LQEND(KQT+NDVMXQ))                     GO TO 998
          IDIV2  = JQDVSY

  510     IF (L2.LT.LQEND(KQT+IDIV2 ))                     GO TO 520
          IDIV2  = IDIV2  + 1
                                                           GO TO 510

  520     IF (L2.LT.LQSTA(KQT+IDIV2 ))                     GO TO 998

#if defined(CERNLIB_QDEVZE)
      IF (NQDEVZ.NE.0) WRITE(IQPRNT,
     +'(''0DEVZE ZSWAP,  Divs : L1,IDIV1,L2,IDIV2= '',4I8)')
     + L1,IDIV1,L2,IDIV2
#endif


          DO 1700 JQDIVI = 1,NDVMXQ
              IF (JQDIVI.GT.JQDVLL.AND.JQDIVI.LT.JQDVSY)   GO TO 1700

              LN = LQSTA(KQT+JQDIVI)
 1300         IF (LN.LT.LQEND(KQT+JQDIVI))       THEN
                  CALL MZCHLN(NCHEKQ,LN)
                  IF (IQFOUL.NE.0)                 THEN
                      WRITE(IQPRNT,
     X                '(''0ZSWAP - BANK INVALID AT'',I10)') LN
                                                           GO TO 998
                  ENDIF
                  LN = IQNX
                  IF (JBIT(IQ(KQS+IQLS),IQDROP).EQ.0)           THEN
                      DO 1400 L = IQLS-IQNL+KQS , IQLS-IQNS+KQS-1
                          IF (LQ(L).EQ.L1) THEN
                              LQ(L) = L2
                          ELSEIF(LQ(L).EQ.L2) THEN
                              LQ(L) = L1
                          ENDIF
 1400                 CONTINUE
                  ENDIF
                                                           GO TO 1300
              ENDIF
 1700     CONTINUE

      ENDIF

                                                           GO TO 999

  998 CONTINUE
#include "zebra/qtofatal.inc"

#include "zebra/qtrace99.inc"
      END
