*
* $Id$
*
* $Log$
* Revision 1.2  1996/04/24 17:26:04  mclareni
* Extend the include file cleanup to dzebra, rz and tq, and also add
* dependencies in some cases.
*
* Revision 1.1.1.1  1996/03/06 10:47:06  mclareni
* Zebra
*
*
*-----------------------------------------------------------
#include "zebra/pilot.h"
#if (defined(CERNLIB_DEBUGON))&&(defined(CERNLIB_VFORT))
#include "zebra/debugvf1.inc"
#endif
      SUBROUTINE DZBKXR(LBANK)
      SAVE NDVXRF
#include "zebra/bankparq.inc"
#include "zebra/divparq.inc"
#include "zebra/questparq.inc"
#include "zebra/mqsys.inc"
#include "zebra/qequ.inc"
#include "zebra/mzcn.inc"
#include "zebra/dzc1.inc"

      CHARACTER CHROUT*(*),CHSTAK*6
      PARAMETER (CHROUT = 'DZBKXR')

#include "zebra/q_jbit.inc"
#include "zebra/q_jbyt.inc"

#if (defined(CERNLIB_DEBUGON))&&(defined(CERNLIB_VFORT))
#include "zebra/debugvf2.inc"
#endif


      CHSTAK          = CQSTAK(MCQSIQ:)
      CQSTAK(MCQSIQ:) = CHROUT

      IF (LBANK.EQ.0) THEN



          IXREF  = IQRTO(KQT+JQDIVI)
          DO 110 I=1,NDVMXQ
              IF (JBIT(IXREF,I).NE.0) THEN
                  IQUEST(30+I) = 1
              ELSE
                  IQUEST(30+I) = 0
              ENDIF
  110     CONTINUE
          IF (JBIT(IXREF,JDVUSQ).EQ.1) THEN
              DO 120 I=1,JQDVLL
  120         IQUEST(30+I) = 1
          ENDIF
          IF (JBYT(IXREF,JDVLGQ,2).NE.0) THEN
              DO 130 I=JQDVSY+1,NDVMXQ
  130         IQUEST(30+I) = 1
          ENDIF
          IF (JBIT(IXREF,JDVSYQ).EQ.1) IQUEST(30+JQDVSY) = 1
          IQUEST(30+JQDIVI) = 1
          NDVXRF = 0
          DO 140 I=1,NDVMXQ
              IF (I.GT.JQDVLL.AND.I.LT.JQDVSY)             GO TO 140
              IF (IQUEST(30+I).EQ.1) THEN
                  NDVXRF = NDVXRF + 1
                  IQUEST(30+NDVXRF) = I
              ENDIF
  140     CONTINUE

                                                           GO TO 999
      ELSE


          LB = LBANK-NL+KQS
          LE = LBANK   +KQS
          DO 220 L=LB,LE
              IF (LQ(L).EQ.0)                              GO TO 220
              DO 210 I=1,NDVXRF
                  IDVXRF = IQUEST(30+I)
                  IF (LQ(L).GE.LQSTA(KQT+IDVXRF).AND.
     X                LQ(L).LT.LQEND(KQT+IDVXRF)    )      GO TO 220
  210         CONTINUE
              CALL DZBKDV(LBANK)
              IF (IQUEST(1).NE.0)                          GO TO 999
              CQINFO = CQDIV//'/'
              CALL DZBKDV(LQ(L))
              IF (IQUEST(1).NE.0)                          GO TO 999
              WRITE(CQINFO(10:),'(I8,'','',A,''/'',I8)')
     X        LBANK,CQDIV,LQ(L)
              CALL DZTEXT(MBKX1Q,CDUMMQ,0)
              IQUEST(1) = 1
                                                           GO TO 999
  220     CONTINUE


          LUP = LQLUP(KQS+LBANK)

          IF (LUP.NE.0) THEN
              DO 310 I=1,NDVXRF
                  IDVXRF = IQUEST(30+I)
                  IF (LUP.GE.LQSTA(KQT+IDVXRF).AND.
     X                LUP.LT.LQEND(KQT+IDVXRF)     )       GO TO 400
  310         CONTINUE
              CALL DZBKDV(LBANK)
              IF (IQUEST(1).NE.0)                          GO TO 999
              CQINFO = CQDIV//'/'
              CALL DZBKDV(LUP)
              IF (IQUEST(1).NE.0)                          GO TO 999
              WRITE(CQINFO(10:),'(I8,'','',A,''/'',I8)')
     X        LBANK,CQDIV,LUP
              CALL DZTEXT(MBKX2Q,CDUMMQ,0)
              IQUEST(1) = 1
                                                           GO TO 999
          ENDIF


  400     LSUP = LQLORG(KQS+LBANK)

          IF (LSUP.EQ.0.AND.LUP.EQ.0)                      GO TO 999
          IF (LSUP.GE.LQSTA(KQT+1).AND.
     X        LSUP.LE.LQSTA(KQT+NDVMXQ+1)) THEN


              DO 410 I=1,NDVXRF
                  IDVXRF = IQUEST(30+I)
                  IF (LSUP.GE.LQSTA(KQT+IDVXRF).AND.
     X                LSUP.LT.LQEND(KQT+IDVXRF)     )      GO TO 999
  410         CONTINUE

              CALL DZBKDV(LBANK)
              IF (IQUEST(1).NE.0)                          GO TO 999
              CQINFO = CQDIV//'/'
              CALL DZBKDV(LSUP)
              IF (IQUEST(1).NE.0)                          GO TO 999
              WRITE(CQINFO(10:),'(I8,'','',A,''/'',I8)')
     X        LBANK,CQDIV,LSUP
              CALL DZTEXT(MBKX3Q,CDUMMQ,0)
              IQUEST(1) = 1
                                                           GO TO 999
          ENDIF
      ENDIF

  999 CQSTAK(MCQSIQ:) = CHSTAK
      RETURN
      END
