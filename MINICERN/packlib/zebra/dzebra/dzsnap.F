*
* $Id$
*
* $Log$
* Revision 1.2  1996/04/24 17:26:15  mclareni
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
      SUBROUTINE DZSNAP (CHTEXT,IXDIV,CHOPT)
#include "zebra/mqsys.inc"
#include "zebra/qequ.inc"
#include "zebra/mzcn.inc"
#include "zebra/zbcdch.inc"
#include "zebra/zbcdk.inc"
#include "zebra/zmach.inc"
#include "zebra/zstate.inc"
#include "zebra/zunit.inc"
#include "zebra/dzc1.inc"
#include "zebra/bankparq.inc"
#include "zebra/divparq.inc"
#include "zebra/questparq.inc"
#include "zebra/storparq.inc"
      CHARACTER *(*) CHTEXT,CHOPT

      CHARACTER CHROUT*(*)
      PARAMETER (CHROUT = 'DZSNAP')

#include "zebra/q_jbit.inc"
#include "zebra/q_jbyt.inc"

#if (defined(CERNLIB_DEBUGON))&&(defined(CERNLIB_VFORT))
#include "zebra/debugvf2.inc"
#endif

      CQSTAK = CHROUT//'/'
      IQUEST(1) = 0
      IF (CHTEXT.NE.CDUMMQ)                                THEN
          CQMAP(1)          = ' '
          CQMAP(2)(1:12)    = ' DZSNAP --- '
          CQMAP(2)(13:100)  = CHTEXT
          CQMAP(2)(101:110) = 'OPTIONS : '
          CQMAP(2)(111:130) = CHOPT
          CALL DZTEXT(0,CDUMMQ,2)
      ENDIF

      CALL DZOPT(CHOPT)


      CALL MZSDIV(IXDIV,0)

      CQMAP(1) = ' '
      CQMAP(2) = '   NAME       LQSTOR NQSTRU  NQREF NQLINK '//
     X           'LQMINR LQ2END JQDVLL JQDVSY NQFEND  LOW-1  LOW-N '//
     X           'HIGH-1 HIGH-N SYST-1 SYST-N    END'
#if defined(CERNLIB_OCTMAP)
      WRITE(CQMAP(3),'(2X,2A4,''('',O8,'')'',15I7)')
#endif
#if !defined(CERNLIB_OCTMAP)
      WRITE(CQMAP(3),'(2X,2A4,''('',Z8,'')'',15I7)')
#endif
     W  NQSNAM(1)     , NQSNAM(2) ,
*  Map  addresses expressed in machine words
#if defined(CERNLIB_WORDMAP)
     W  LQSTOR+1 ,
*  Map addresses expressed in bytes
#endif
#if !defined(CERNLIB_WORDMAP)
     W  (LQSTOR+1)*4 ,
#endif
     W  NQSTRU , NQREF    , NQLINK , NQMINR        ,
     W  LQ2END , JQDVLL   , JQDVSY , NQFEND        ,
     W  LQSTA (KQT+MDVLWQ),LQEND (KQT+MDVLWQ)-1    ,
     W  LQSTA (KQT+MDVHGQ),LQEND (KQT+MDVHGQ)-1    ,
     W  LQSTA (KQT+JQDVSY),LQEND (KQT+JQDVSY)-1    ,
     W  LQSTA (KQT+NDVMXQ+1)-1

      CALL DZTEXT(0,CDUMMQ,3)



      DO 100 IFENCE=-NQFEND+1,0
          IF(LQ(KQS+IFENCE).NE.IQNIL)      THEN
              WRITE (CQINFO,'(I5,1X,Z16)') IFENCE,LQ(KQS+IFENCE)
              CALL DZTEXT(MSNA1Q,CDUMMQ,0)
              IF (IFLOPT(MPOSSQ).NE.0)                     GO TO 998
          ENDIF
  100 CONTINUE


      IF ((LQ(KQS+LQSTA(KQT+21)+1).NE.IQNIL) .OR.
     X    (LQ(KQS+LQSTA(KQT+21)+2).NE.IQNIL)      ) THEN
          WRITE (CQINFO,'(Z16,1X,Z16)')
     X    LQ(KQS+LQSTA(KQT+21)+1),LQ(KQS+LQSTA(KQT+21)+2)
          CALL DZTEXT(MSNA2Q,CDUMMQ,0)
          IF (IFLOPT(MPOSSQ).NE.0)                         GO TO 998
      ENDIF

      IF(NQSTRU.GT.NQREF)              THEN
          WRITE (CQINFO,'(I8,''>'',I8)') NQSTRU,NQREF
          CALL DZTEXT(MSNA3Q,CDUMMQ,0)
          IF (IFLOPT(MPOSSQ).NE.0)                         GO TO 998
      ENDIF

      IF(NQREF.GT.NQLINK)              THEN
          WRITE (CQINFO,'(I8,''>'',I8)') NQREF,NQLINK
          CALL DZTEXT(MSNA4Q,CDUMMQ,0)
          IF (IFLOPT(MPOSSQ).NE.0)                         GO TO 998
      ENDIF

      IF(LQSTA(KQT+2)-LQEND(KQT+1).LT.NQMINR)      THEN
          WRITE (CQINFO,'(I8,''-'',I8,''<'',I8)')
     X     LQSTA(KQT+2),LQEND(KQT+1),NQMINR
          CALL DZTEXT(MSNA5Q,CDUMMQ,0)
          IF (IFLOPT(MPOSSQ).NE.0)                         GO TO 998
      ENDIF

      IF(NQMINR.GT.LQ2END)             THEN
          WRITE (CQINFO,'(I8,''>'',I8)') NQMINR,LQ2END
          CALL DZTEXT(MSNA6Q,CDUMMQ,0)
          IF (IFLOPT(MPOSSQ).NE.0)                         GO TO 998
      ENDIF

      IF(LQ2END.GT.LQSTA(KQT+21))      THEN
          WRITE (CQINFO,'(I8,''>'',I8)') LQ2END,LQSTA(KQT+21)
          CALL DZTEXT(MSNA7Q,CDUMMQ,0)
          IF (IFLOPT(MPOSSQ).NE.0)                         GO TO 998
      ENDIF


      IF (IFLOPT(MPOSTQ).NE.0)         THEN
          IFLOPT(MPOSCQ) = 1
          CALL UCOPY(IFLOPT,IQUEST(71),26)
          CALL VZERO(IFLOPT,26)
          IFLOPT(MPOSNQ) = 1
          IFLOPT(MPOSQQ) = 1
          IFLOPT(MPOSTQ) = 1
          CALL DZARE1('DZSNAP  L option',' ',0,'NQT')
          CALL UCOPY(IQUEST(71),IFLOPT,26)
      ENDIF
      IF (IFLOPT(MPOSLQ).NE.0) THEN
          CALL UCOPY(IFLOPT,IQUEST(71),26)
          CALL VZERO(IFLOPT,26)
          IFLOPT(MPOSNQ) = 1
          CALL DZARE1('DZSNAP  L option',' ',0,'N')
          CALL UCOPY(IQUEST(71),IFLOPT,26)
      ENDIF
      IF (IFLOPT(MPOSWQ).NE.0)               THEN
          CQMAP(1) = ' '
#if !defined(CERNLIB_OCTMAP)
          WRITE(CQMAP(2),'(''  WORKING SPACE   ADR(LQ(0)) = '',Z8)')
#endif
#if defined(CERNLIB_OCTMAP)
          WRITE(CQMAP(2),'(''  WORKING SPACE   ADR(LQ(0)) = '',O8)')
*  Map  addresses expressed in machine words
#endif
#if defined(CERNLIB_WORDMAP)
     W      LQSTOR+1
*  Map addresses expressed in bytes
#endif
#if !defined(CERNLIB_WORDMAP)
     W      (LQSTOR+1)*4
#endif
          CALL DZTEXT(0,CDUMMQ,2)
          LBASE  = 0
          IBASE  = 0
          JDFD   = NQLINK + 1
          NDW    = LQSTA(KQT+1) - 1
          IF (IFLOPT(MPOSTQ).NE.0) THEN
              NDW = MIN(NDW,NQWCUT+NQLINK)
          ELSE
              NDW = MIN(NDW,NQLINK)
          ENDIF
          CALL DZDATA(CDUMMQ)
      ENDIF


      IF (JBYT(IXDIV,1,JSTIDQ-1).EQ.0)  THEN
          JJDIV = MZIXCO(IXDIV+21,IXDIV+22,0,0)
          JJDIV = MZDVAC(JJDIV)
      ELSE
          JJDIV  = MZDVAC (IXDIV)
      ENDIF

      IF (IFLOPT(MPOSEQ)+IFLOPT(MPOSFQ)+IFLOPT(MPOSMQ).EQ.0) GO TO 999



      NDZRSV = 0
      CALL DZBKUP(0)
      IF (IQUEST(1).NE.0)                                  GO TO 999

      DO 1000 JDIVI = 1,NDVMXQ
          IF ( JDIVI.GT.JQDVLL.AND.JDIVI.LT.JQDVSY)        GO TO 1000
          IF (JBIT(JJDIV,JDIVI).EQ.0)                      GO TO 1000
          WRITE(CQMAP,'(1X,/,
     W       '' DZSNAP.   -----  Store nb.'',I2,'' = '',2A4,
     W       '' Division nb.'',I2,'' = '',2A4,20X,20(''-''),/)')
     W      JQSTOR,NQSNAM(1),  NQSNAM(2),
     W      JDIVI,IQDN1(KQT+JDIVI),IQDN2(KQT+JDIVI)
          CALL DZTEXT(0,CDUMMQ,3)
          LN     = LQSTA(KQT+JDIVI)
          LSTOP  = LQEND(KQT+JDIVI)
          IF(LN.EQ.LSTOP)                  THEN
              CQLINE = '          -- Division contains no banks --'
              CALL DZTEXT(0,CDUMMQ,1)
                                                           GO TO 1000
          ENDIF



          JQDIVI = JDIVI
          CALL DZBKXR(0)
          IF (IQUEST(1).NE.0)                              GO TO 999

  300     IF(LN.LT.LSTOP)                  THEN
              CALL DZMAP
              IF (IQUEST(1).NE.0)                          GO TO 400
              LN  = LX
                                                           GO TO 300
          ENDIF

                                                           GO TO 1000
  400     WRITE(CQINFO,'(I8)') LN
          CALL DZTEXT(MSNA8Q,CDUMMQ,0)
          LBKCL= LN
          IFLOPT(MPOSWQ) = -1
  500     LBK = LN
  600     LBK = LBK + 1
          IF (LBK.GE.LSTOP)                THEN
              LN     = LSTOP
              LBASE  = LBKCL-1
              IBASE  = LBASE
              NDW    = MIN(LN - LBKCL,NQWCUT)
              JDFD   = 1
              CALL DZDATA(CDUMMQ)
                                                           GO TO 1000
          ENDIF
          CALL MZCHLN (NCHEKQ,LBK)
          IF (IQFOUL.NE.0)                                 GO TO 600
          LN     = LBK
          LBK    = IQNX
          IF (LBK.GE.LSTOP)                THEN
              LN     = LSTOP
              LBASE  = LBKCL-1
              IBASE  = LBASE
              NDW    = MIN(LN - LBKCL,NQWCUT)
              JDFD   = 1
              CALL DZDATA(CDUMMQ)
                                                           GO TO 1000
          ENDIF
          CALL MZCHLN (NCHEKQ,LBK)
          IF (IQFOUL.NE.0)                                 GO TO 500
  700     CONTINUE
              LBASE  = LBKCL-1
              IBASE  = LBASE
              NDW    = MIN(LN - LBKCL,NQWCUT)
              JDFD   = 1
              CALL DZDATA(CDUMMQ)
          WRITE(CQMAP,'(1X,/,'' RECOVER AT ADR'',I8)')   LN
          CALL DZTEXT(0,CDUMMQ,2)
          IQUEST(1) = 0
                                                           GO TO 300

 1000 CONTINUE

                                                           GO TO 999

  998 IQUEST(1) = 1

  999 RETURN
      END
