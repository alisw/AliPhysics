*
* $Id$
*
* $Log$
* Revision 1.2  1996/04/24 17:26:20  mclareni
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
      SUBROUTINE DZVERI (CHTEXT,IXDIV,CHOPT)
#include "zebra/mqsys.inc"
#include "zebra/qequ.inc"
#include "zebra/mzcn.inc"
#include "zebra/zbcdch.inc"
#include "zebra/zbcdk.inc"
#include "zebra/zmach.inc"
#include "zebra/zunit.inc"
#include "zebra/dzc1.inc"
#include "zebra/bankparq.inc"
#include "zebra/divparq.inc"
#include "zebra/questparq.inc"
#include "zebra/storparq.inc"

      PARAMETER (KINUSQ=2**(IDVUSQ-1),KINSYQ=2**(IDVSYQ-1))

      CHARACTER*(*) CHOPT,CHTEXT
      CHARACTER CHROUT*(*)
      PARAMETER (CHROUT = 'DZVERI')

#include "zebra/q_jbit.inc"
#include "zebra/q_jbyt.inc"

#if (defined(CERNLIB_DEBUGON))&&(defined(CERNLIB_VFORT))
#include "zebra/debugvf2.inc"
#endif

      CQSTAK = CHROUT//'/'
      IQUEST(1) = 0

      CALL MZSDIV (IXDIV,-1)

      CALL DZOPT(CHOPT)



      ICHAIN = 1
      ISTRUC = 0

      IF (IFLOPT(MPOSSQ).NE.0.AND.IFLOPT(MPOSCQ).EQ.0) THEN
          ICHAIN = 0
      ENDIF
      IF (IFLOPT(MPOSLQ).NE.0)          THEN
          ICHAIN = 1
          ISTRUC = 1
      ENDIF

      IF (IFLOPT(MPOSSQ).EQ.0)                             GO TO 1000


      DO 100 IFENCE=-NQFEND+1,0
          IF(LQ(KQS+IFENCE).NE.IQNIL)      THEN
              WRITE (CQINFO,'(I5,1X,Z16)') IFENCE,LQ(KQS+IFENCE)
              CALL DZTEXT(MSNA1Q,CDUMMQ,0)
                                                           GO TO 998
          ENDIF
  100 CONTINUE


      IF ((LQ(KQS+LQSTA(KQT+21)+1).NE.IQNIL) .OR.
     X    (LQ(KQS+LQSTA(KQT+21)+2).NE.IQNIL)      ) THEN
          WRITE (CQINFO,'(Z16,1X,Z16)')
     X    LQ(KQS+LQSTA(KQT+21)+1),LQ(KQS+LQSTA(KQT+21)+2)
          CALL DZTEXT(MSNA2Q,CDUMMQ,0)
                                                           GO TO 998
      ENDIF

      IF(NQSTRU.GT.NQREF)              THEN
          WRITE (CQINFO,'(I8,''>'',I8)') NQSTRU,NQREF
          CALL DZTEXT(MSNA3Q,CDUMMQ,0)
                                                           GO TO 998
      ENDIF

      IF(NQREF.GT.NQLINK)              THEN
          WRITE (CQINFO,'(I8,''>'',I8)') NQREF,NQLINK
          CALL DZTEXT(MSNA4Q,CDUMMQ,0)
                                                           GO TO 998
      ENDIF

      IF(LQSTA(KQT+2)-LQEND(KQT+1).LT.NQMINR)      THEN
          WRITE (CQINFO,'(I8,''-'',I8,''<'',I8)')
     X     LQSTA(KQT+2),LQEND(KQT+1),NQMINR
          CALL DZTEXT(MSNA5Q,CDUMMQ,0)
                                                           GO TO 998
      ENDIF

      IF(NQMINR.GT.LQ2END)             THEN
          WRITE (CQINFO,'(I8,''>'',I8)') NQMINR,LQ2END
          CALL DZTEXT(MSNA6Q,CDUMMQ,0)
                                                           GO TO 998
      ENDIF

      IF(LQ2END.GT.LQSTA(KQT+21))      THEN
          WRITE (CQINFO,'(I8,''>'',I8)') LQ2END,LQSTA(KQT+21)
          CALL DZTEXT(MSNA7Q,CDUMMQ,0)
                                                           GO TO 998
      ENDIF



      LPREV = NQLINK
      DO 300 JDIVI=1,NDVMXQ
          IF(JDIVI.GT.JQDVLL.AND.JDIVI.LT.JQDVSY)          GO TO 300

          CALL UHTOC (IQDN1(KQT+JDIVI),4,CQDIV(1:4),4)
          CALL UHTOC (IQDN2(KQT+JDIVI),4,CQDIV(5:8),4)
          CQINFO = CQDIV//'/'

          IF(LQSTA(KQT+JDIVI).GT.LQEND(KQT+JDIVI))       THEN
              WRITE (CQINFO(10:),'(I8,''>'',I8)')
     X         LQSTA(KQT+JDIVI),LQEND(KQT+JDIVI)
              CALL DZTEXT(MVER1Q,CDUMMQ,0)
                                                           GO TO 998
          ENDIF

          IF(LQSTA(KQT+JDIVI).LT.LPREV            )       THEN
              WRITE (CQINFO(10:),'(I8,''<'',I8)') LQSTA(KQT+JDIVI),LPREV
              CALL DZTEXT(MVER2Q,CDUMMQ,0)
                                                           GO TO 998
          ENDIF

          IF((LQEND(KQT+JDIVI)-LQSTA(KQT+JDIVI)).GT.
     X        NQDMAX(KQT+JDIVI)                       )    THEN
              WRITE (CQINFO(10:),'(I8,''>'',I8)')
     X        (LQEND(KQT+JDIVI)-LQSTA(KQT+JDIVI)),NQDMAX(KQT+JDIVI)
              CALL DZTEXT(MVER3Q,CDUMMQ,0)
                                                           GO TO 998
          ENDIF

          IF(IQMODE(KQT+JDIVI).NE.IDVFWQ .AND.
     X       IQMODE(KQT+JDIVI).NE.IDVBWQ     )            THEN
              WRITE (CQINFO(10:),'(I10)') IQMODE(KQT+JDIVI)
              CALL DZTEXT(MVER4Q,CDUMMQ,0)
                                                           GO TO 998
          ENDIF

          MKIND = IQKIND(KQT+JDIVI)

          DO 200 I=1,NDVMXQ
              IF(I.NE.JDIVI.AND.JBIT(MKIND,I).EQ.1)       THEN
              WRITE (CQINFO(10:),'(I2,'',        /'',I2)') JDIVI,I
              CALL UHTOC (IQDN1(KQT+JDIVI),4,CQINFO(13:16),4)
              CALL UHTOC (IQDN2(KQT+JDIVI),4,CQINFO(17:20),4)
              CALL DZTEXT(MVER5Q,CDUMMQ,0)
                                                           GO TO 998
          ENDIF
  200     CONTINUE

          MKIND = JBYT(MKIND,JDVUSQ,JDVSYQ-NDVIDQ)
          IF(JDIVI.LE.JQDVLL)                               THEN
              IF(MKIND.NE.KINUSQ)                              THEN
              WRITE (CQINFO(10:),'(4I5)') JDIVI,JQDVLL,JQDVSY,MKIND
              CALL DZTEXT(MVER6Q,CDUMMQ,0)
                                                           GO TO 998
              ENDIF
          ELSEIF(JDIVI.EQ.JQDVSY)                           THEN
              IF(MKIND.NE.KINSYQ)                              THEN
              WRITE (CQINFO(10:),'(4I5)') JDIVI,JQDVLL,JQDVSY,MKIND
              CALL DZTEXT(MVER7Q,CDUMMQ,0)
                                                           GO TO 998
              ENDIF
          ELSE
              IF(MKIND.EQ.KINUSQ.OR.MKIND.EQ.KINSYQ)           THEN
              WRITE (CQINFO(10:),'(4I5)') JDIVI,JQDVLL,JQDVSY,MKIND
              CALL DZTEXT(MVER8Q,CDUMMQ,0)
                                                           GO TO 998
              ENDIF
          ENDIF
          LPREV = LQEND(KQT+JDIVI)
  300 CONTINUE



 1000 IF (ICHAIN.EQ.0)                                     GO TO 1999

      IF (ISTRUC.NE.0) THEN


          NDZRSV = 0
          CALL DZBKUP(0)
          IF (IQUEST(1).NE.0)                              GO TO 998

      ENDIF


      CALL UCOPY(IFLOPT,IQUEST(71),26)
      CALL VZERO(IFLOPT,26)
      IFLOPT(MPOSNQ) = 1
      IFLOPT(MPOSQQ) = 1
      CALL DZARE1('DZVERI  L option',' ',0,'NQ')
      CALL UCOPY(IQUEST(71),IFLOPT,26)


      IF (JBYT(IXDIV,1,JSTIDQ-1).EQ.0)  THEN
          JJDIV = MZIXCO(IXDIV+21,IXDIV+22,0,0)
          JJDIV = MZDVAC(JJDIV)
      ELSE
          JJDIV  = MZDVAC (IXDIV)
      ENDIF

      DO 1700 JDIVI = 1,NDVMXQ
          IF (JDIVI.GT.JQDVLL.AND.JDIVI.LT.JQDVSY)         GO TO 1700

          IF (JBIT(JJDIV,JDIVI).EQ.0)                      GO TO 1700

          CALL UHTOC (IQDN1(KQT+JDIVI),4,CQDIV(1:4),4)
          CALL UHTOC (IQDN2(KQT+JDIVI),4,CQDIV(5:8),4)
          CQINFO = CQDIV//'/'


          JQDIVI = JDIVI
          CALL DZBKXR(0)
          IF (IQUEST(1).NE.0)                              GO TO 998


          LN = LQSTA(KQT+JDIVI)
 1300     IF (LN.LT.LQEND(KQT+JDIVI))       THEN
              CALL MZCHLN(NCHEKQ,LN)
              IF (IQFOUL.NE.0)   THEN
                  WRITE(CQINFO(10:),'(3I10)') JDIVI,LN,IQFOUL
                  CALL DZTEXT(MVER9Q,CDUMMQ,0)
                                                           GO TO 998
              ENDIF
              LN = IQNX
              LS = IQLS
              NL = IQNL
              IF (IQND.LT.0)                               GO TO 1300
              IF (JBIT(IQ(KQS+LS),IQDROP).EQ.1)            GO TO 1300


              IF (ISTRUC.EQ.1)                 THEN
                  LB = IQLS-IQNS+KQS
                  LE = IQLS     +KQS
                  DO 1400 L=LB,LE
                      IF (LQ(L).EQ.0)                      GO TO 1400
                      CALL MZCHLS(NCHEKQ,LQ(L))
                      IF (IQFOUL.NE.0) THEN
                          CALL DZBKDV(LQ(L))
                          IF (IQUEST(1).NE.0)              GO TO 998
                          WRITE(CQINFO(10:),
     X                      '(I8,'','',A,''/'',I8,'','',I4)')
     X                       L-KQS,CQDIV,LQ(L),IQFOUL
                                                           GO TO 998
                      ENDIF
 1400             CONTINUE
                  CALL DZBKUP(LS)
                  IF (IQUEST(1).NE.0)                      GO TO 998
                  CALL DZBKXR(LS)
                  IF (IQUEST(1).NE.0)                      GO TO 998
              ENDIF

                                                           GO TO 1300
          ENDIF
 1700 CONTINUE

 1999 IF (CHTEXT.NE.CDUMMQ)                 THEN
          CQMAP(1)          = ' '
          CQMAP(2)          = ' DZVERI --- '
          CQMAP(2)(13:99)   = CHTEXT
          CQMAP(2)(101:110) = 'OPTIONS : '
          CQMAP(2)(111:120) = CHOPT
          CQMAP(2)(126:)    = 'OK'
          CALL DZTEXT(0,CDUMMQ,2)
      ENDIF
                                                           GO TO 999


  998 IF (CHTEXT.NE.CDUMMQ)                 THEN
          CQMAP(1)          = ' '
          CQMAP(2)          = ' DZVERI --- '
          CQMAP(2)(13:99)   = CHTEXT
          CQMAP(2)(101:110) = 'OPTIONS : '
          CQMAP(2)(111:115) = CHOPT
          CQMAP(2)(117:)    = '??PROBLEMS ? '
          CALL DZTEXT(0,CDUMMQ,2)
      ENDIF
      IQUEST(1)=1


      NQFATA = 2
      IQUEST(11) = JQSTOR
      IQUEST(12) = JDIVI
      IF(ICHAIN.NE.0)      THEN
          IQUEST(13) = LN
          IQUEST(14) = IQLS
          IQUEST(15) = IQNL
          IQUEST(16) = IQNS
          IQUEST(17) = IQND
          NQFATA = 17
          IF(ISTRUC.EQ.1)                            THEN
              NQFATA = 21
              IF(L.GE.LQSTA(KQT+JDIVI).AND.L.LT.LQEND(KQT+JDIVI)) THEN
                  IQUEST(18) = L
                  IQUEST(19) = LQ(L)
              ELSE
                  IQUEST(18) = 0
                  IQUEST(19) = 0
              ENDIF
              IQUEST(20) = LQLUP(KQS+LS)
              IQUEST(21) = LQLORG(KQS+LS)
          ENDIF
      ENDIF


      IF(IFLOPT(MPOSFQ).NE.0) CALL ZFATAM(CHROUT)


  999 RETURN
      END
