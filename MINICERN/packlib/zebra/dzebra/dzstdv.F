*
* $Id$
*
* $Log$
* Revision 1.2  1996/04/24 17:26:16  mclareni
* Extend the include file cleanup to dzebra, rz and tq, and also add
* dependencies in some cases.
*
* Revision 1.1.1.1  1996/03/06 10:47:07  mclareni
* Zebra
*
*
*-----------------------------------------------------------
#include "zebra/pilot.h"
#if (defined(CERNLIB_DEBUGON))&&(defined(CERNLIB_VFORT))
#include "zebra/debugvf1.inc"
#endif
      SUBROUTINE DZSTDV
      SAVE CKINDV,CMODV
#include "zebra/bankparq.inc"
#include "zebra/divparq.inc"
#include "zebra/mqsys.inc"
#include "zebra/qequ.inc"
#include "zebra/mzcn.inc"
#include "zebra/zunit.inc"
#include "zebra/dzc1.inc"
      CHARACTER CKINDV(4)*7,CMODV(0:1)*5

      CHARACTER CHROUT*(*),CHSTAK*6
      PARAMETER (CHROUT = 'DZSTDV')

#include "zebra/q_jbit.inc"
#include "zebra/q_jbyt.inc"

#if (defined(CERNLIB_DEBUGON))&&(defined(CERNLIB_VFORT))
#include "zebra/debugvf2.inc"
#endif

      DATA CKINDV / 'U/EVENT' , ' U/LONG' , 'PACKAGE' , ' SYSTEM' /
      DATA CMODV / 'FORWD'   , 'REVRS'   /


      CHSTAK          = CQSTAK(MCQSIQ:)
      CQSTAK(MCQSIQ:) = CHROUT


      CQMAP(1) = ' '
      CQMAP(2) = '  --- Division parameters ---'
      CQMAP(3) = ' '
      CQMAP(4) = '    DIVISION    START    END       MAX    '//
     X           'KIND   MODE  WIPES  GARB.  GARB. PUSHES'//
     X           '      LIVE BANKS  DROPPED BANKS    BANKS TOTAL'
      CQMAP(5) = '  NB.   NAME   ADDRESS ADDRESS  LENGTH    '//
     X           '                    SYST.   FREE       '//
     X           '  NUMB.   LENGTH NUMB.   LENGTH NUMB.   LENGTH'
      CQMAP(6) = ' ========================================='//
     X           '======================================='//
     X           '=============================================='
      CALL DZTEXT(0,CDUMMQ,6)


      DO 300 JDIVI = 1,NDVMXQ


          IF (JDIVI.GT.JQDVLL.AND.JDIVI.LT.JQDVSY)         GO TO 300
          IMODE  = JBIT(IQMODE(KQT+JDIVI),JDVBFQ)
          JKINDU = JBYT(IQKIND(KQT+JDIVI),JDVUSQ,JDVSYQ-NDVIDQ)
          IKIND  = MIN(JKINDU/2+1,IDVSYQ)

             NBANK  = 0
             NDROP  = 0
             NBANKL = 0
             NDROPL = 0
             LN     = LQSTA (KQT+JDIVI)
             LIMIT  = LQEND (KQT+JDIVI)
             IQFOUL = 0

100          IF (LN.GE.LIMIT)                              GO TO 200
                 CALL MZCHLN(NCHEKQ,LN)
                 IF(IQFOUL.NE.0)                           GO TO 200
                 NBKTOT = IQNX - LN

                 IF     (IQND.LT.0)                          THEN
                     NDROP  = NDROP  + 1
                     NDROPL = NDROPL + NBKTOT
                 ELSE
                     IF     (JBIT(IQ(KQS+IQLS),IQDROP).EQ.1)     THEN
                         NDROP  = NDROP  + 1
                         NDROPL = NDROPL + NBKTOT
                     ELSE
                         NBANK  = NBANK  + 1
                         NBANKL = NBANKL + NBKTOT
                     ENDIF
                 ENDIF

                 LN     = IQNX

                                                           GO TO 100

  200     IF (IQFOUL.EQ.0)                            THEN
              WRITE(CQLINE,
#if !defined(CERNLIB_QMIBM)
     W        '(I4,2X,2A4,3I8,1X,A7,2X,A5,4I7,2X,3(I6,I9))',ERR=300)
#endif
#if defined(CERNLIB_QMIBM)
     W        '(I4,2X,2A4,3I8,1X,A7,2X,A5,4I7,2X,3(I6,I9))')
#endif
     W        JDIVI , IQDN1(KQT+JDIVI) , IQDN2(KQT+JDIVI) ,
     W        LQSTA(KQT+JDIVI),LQEND(KQT+JDIVI)-1,NQDSIZ(KQT+JDIVI),
     W        CKINDV(IKIND) , CMODV(IMODE) , NQDWIP(KQT+JDIVI) ,
     W        NQDGAU(KQT+JDIVI),NQDGAF(KQT+JDIVI),NQDPSH(KQT+JDIVI),
     W        NBANK , NBANKL , NDROP , NDROPL ,
     W        NBANK + NDROP  , NBANKL + NDROPL
              CALL DZTEXT(0,CDUMMQ,1)

          ELSE


              WRITE(CQMAP,'('' ?????? Bank beginning at address='',I8,
     W          '' has problems   -- no bank statistics'',/,
#if !defined(CERNLIB_QMIBM)
     W          I4,2X,2A4,3I8,1X,A7,2X,A5,4I7,2X)',ERR=300)
#endif
#if defined(CERNLIB_QMIBM)
     W          I4,2X,2A4,3I8,1X,A7,2X,A5,4I7,2X)')
#endif
     W        LN,           JDIVI,IQDN1(KQT+JDIVI),IQDN2(KQT+JDIVI),
     W        LQSTA(KQT+JDIVI),LQEND(KQT+JDIVI)-1,NQDSIZ(KQT+JDIVI),
     W        CKINDV(IKIND) , CMODV(IMODE) , NQDWIP(KQT+JDIVI) ,
     W        NQDGAU(KQT+JDIVI),NQDGAF(KQT+JDIVI),NQDPSH(KQT+JDIVI)
              CALL DZTEXT(0,CDUMMQ,2)

          ENDIF

  300 CONTINUE

  999 CQSTAK(MCQSIQ:) = CHSTAK
      END
