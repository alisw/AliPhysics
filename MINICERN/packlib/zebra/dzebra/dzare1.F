*
* $Id$
*
* $Log$
* Revision 1.2  1996/04/24 17:26:00  mclareni
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
      SUBROUTINE DZARE1  (CHTEXT,CLA,LLA,CHOPT)
      SAVE CLATYP
#include "zebra/mqsys.inc"
#include "zebra/qequ.inc"
#include "zebra/mzcn.inc"
#include "zebra/zbcdch.inc"
#include "zebra/zbcdk.inc"
#include "zebra/zstate.inc"
#include "zebra/zunit.inc"
#include "zebra/dzc1.inc"
#include "zebra/questparq.inc"
#include "zebra/storparq.inc"
      CHARACTER *(*) CLA,CHOPT,CHTEXT
      INTEGER ILANAM(2)

      CHARACTER CHROUT*(*),CHSTAK*6,   CLATYP(0:1)*9, CLA8*8, CAKTIV*8
      PARAMETER (CHROUT = 'DZARE1')
      DATA      CLATYP               /'PERMANENT','TEMPORARY'/

#include "zebra/q_jbit.inc"
#include "zebra/q_jbyt.inc"

#if (defined(CERNLIB_DEBUGON))&&(defined(CERNLIB_VFORT))
#include "zebra/debugvf2.inc"
#endif


      CHSTAK          = CQSTAK(MCQSIQ:)
      CQSTAK(MCQSIQ:) = CHROUT

      IF (LQSTA(KQT+JQDVSY).EQ.LQEND(KQT+JQDVSY)) THEN
          IQUEST(1) = 0
                                                           GO TO 999
      ENDIF

      LSYSB  = LQSYSS(KQT+MSYLAQ)
      CALL MZCHLS(NCHEKQ,LSYSB)
      IF (IQFOUL.NE.0)  THEN
          WRITE(CQINFO,'(''Store = '',2A4)') NQSNAM(1),NQSNAM(2)
          CALL DZTEXT(MARE1Q,CHTEXT,0)
                                                           GO TO 999
      ENDIF


      NWTAB  = IQ(KQS+LSYSB+MLAUSQ)
      IF(NWTAB.LE.1)                   THEN
          WRITE(CQINFO,'(''Store = '',2A4)') NQSNAM(1),NQSNAM(2)
          CALL DZTEXT(MARE2Q,CHTEXT,0)
                                                           GO TO 999
      ENDIF

      IF (IFLOPT(MPOSNQ).NE.0)  THEN


          CLA8 = CLA
          CALL UCTOH (CLA8,ILANAM,4,8)
      ELSE


          LLINK = LOCF(LLA) - LQSTOR
      ENDIF


      LENTRY = LSYSB + KQS + MLAUSQ

      IFOUND = 0

      DO 100 IENTRY = 1,(NWTAB-1)/NLAENQ


          IF (IFLOPT(MPOSNQ).NE.0) THEN
              IF(CLA.NE.' ')   THEN
                  IF (ILANAM(1).NE.IQ(LENTRY+MLAN1Q) .OR.
     X                ILANAM(2).NE.IQ(LENTRY+MLAN2Q)  )    GO TO 100
              ELSE
                  IF(IENTRY.LE.2)                          GO TO 100
              ENDIF

          ELSE


              IF (LLINK.LT.IQ(LENTRY+MLAADQ) .OR.
     X            LLINK.GT.IQ(LENTRY+MLALTQ)       )       GO TO 100

          ENDIF


          LLAAR1 = IQ(LENTRY+MLAADQ)
          LLAARL = IQ(LENTRY+MLALTQ)
          NLANS  = IQ(LENTRY+MLANSQ)
          JTEMP  = JBIT(NLANS,JLATMQ)
          NTEMP  = NLATMQ*JTEMP
          NLANS  = JBYT(NLANS,JLANSQ,NLANSQ)-NTEMP
          IF(JTEMP.EQ.0.OR.
     +     (JTEMP.EQ.1.AND.LQ(KQS+LLAAR1+MLACTQ-1).NE.0))          THEN
              CAKTIV = '  ACTIVE'
          ELSE
              CAKTIV = 'INACTIVE'
          ENDIF


          DO 50 I = NTEMP,NLANS-1
              LS  = LQ(KQS+LLAAR1+I)
              IF (LS.EQ.0)                                 GO TO 50
              CALL MZCHLS(NCHEKQ,LS)
              ID     = IQID
              IF (JBIT(IQ(KQS+LS),IQDROP).EQ.1)            GO TO 50
              IF (IQFOUL.EQ.0) THEN
                  LUP    = LQLUP(KQS+LS)
                  IF (LUP.EQ.0)                            GO TO 40
                  CALL MZCHLS(NCHEKQ,LUP)
                  IF (IQFOUL.NE.0) THEN
                      WRITE(CQINFO,
     X                 '(2A4,''/'',I5,''= '',A4,2I10)')
     X                 IQ(LENTRY+MLAN1Q),IQ(LENTRY+MLAN2Q),
     X                 I+1,ID,LS,LUP
                      CALL DZTEXT(MARE3Q,CDUMMQ,0)
                  ENDIF
   40             LORIG  = LQLORG(KQS+LS)
                  IF (LORIG.EQ.0)                          GO TO 50
                  IF(LORIG.LT.IQTABV(KQT+13).OR.LORIG.GT.IQTABV(KQT+14))
     X                  THEN
                      WRITE(CQINFO,'(2A4,''/'',I5,''= '',A4,2I10)')
     X                 IQ(LENTRY+MLAN1Q),IQ(LENTRY+MLAN2Q),
     X                 I+1,ID,LS,LORIG
                      CALL DZTEXT(MARE4Q,CDUMMQ,0)
                  ELSEIF (LQ(KQS+LORIG).NE.LS) THEN
                      WRITE(CQINFO,'(2A4,''/'',I5,''= '',A4,2I10)')
     X                 IQ(LENTRY+MLAN1Q),IQ(LENTRY+MLAN2Q),
     X                 I+1,ID,LS,LORIG
                      CALL DZTEXT(MARE4Q,CDUMMQ,0)
                  ENDIF
                  IF (IQND.LT.NQDCUT.AND.IFLOPT(MPOSTQ).NE.0)
     X             CALL SBIT1 (IQ(KQS+LS),IQCRIT)
              ELSE
                  WRITE(CQINFO,'(2A4,''/'',I5,''='',I10)')
     X             IQ(LENTRY+MLAN1Q),IQ(LENTRY+MLAN2Q),I+1,LS
                  CALL DZTEXT(MARE5Q,CDUMMQ,0)
              ENDIF
   50     CONTINUE

          IF (IFLOPT(MPOSQQ).EQ.0) THEN
             IF (CHTEXT.NE.CDUMMQ)                                 THEN
                 CQMAP(1)      = ' '
                 CQMAP(2)      = ' DZAREA -- '
                 CQMAP(2)(12:) = CHTEXT
                 WRITE(CQMAP(2)(80:),'('' -- Dump of link area '',2A4,
     W           '' Options: '',A)')
     W           IQ(LENTRY+MLAN1Q),IQ(LENTRY+MLAN2Q),CHOPT
                 CALL DZTEXT(0,CDUMMQ,2)
             ENDIF

             CQMAP(1) = ' '
             WRITE(CQMAP(2),
#if !defined(CERNLIB_OCTMAP)
     W       '('' This '',A9,'' LINK AREA is at absolute address '',Z8,
#endif
#if defined(CERNLIB_OCTMAP)
     W       '('' This '',A9,'' LINK AREA is at absolute address '',O11,
#endif
     W       '' with NL/NS ='',I7,''/'',I7,4X,'' and is '',A8)'      )
     W        CLATYP(JTEMP),(LLAAR1+LQSTOR)
#if !defined(CERNLIB_WORDMAP)
     W       *4
#endif
     W        ,LLAARL-LLAAR1-NTEMP,NLANS,CAKTIV
             CALL DZTEXT(0,CDUMMQ,2)
             CALL DZTEXT(1,CDUMMQ,1)


             LBASE  = LLAAR1 + NTEMP - 1
             IBASE  = 0
             NDW    = LLAARL - LLAAR1  - NTEMP
             JDFD   = NDW    + 1

             CALL DZDATA(CDUMMQ)

             IFOUND = 1
             IF (IFLOPT(MPOSNQ).NE.0) THEN
                 IF (CLA.NE.' ')                           GO TO 999
             ELSE
                                                           GO TO 999
             ENDIF
          ELSE
              IQUEST(10) = IENTRY
              IQUEST(11) = IQ(LENTRY+MLAN1Q)
              IQUEST(12) = IQ(LENTRY+MLAN2Q)
              IQUEST(13) = LLAAR1
              IQUEST(14) = LLAAR1 + NLANS
              IQUEST(15) = LLAARL
              IQUEST(16) = NLANS
              IQUEST(17) = LLAARL-LLAAR1
              IQUEST(18) = JTEMP
              IQUEST(19) = LQ(KQS+LLAAR1+MLACTQ-1)
                                                           GO TO 999
             ENDIF
  100     LENTRY = LENTRY + NLAENQ

      IF (IFOUND.EQ.0)           THEN
          IF (IFLOPT(MPOSQQ).EQ.0) THEN
             IF (IFLOPT(MPOSNQ).NE.0) THEN
                 WRITE(CQINFO,'(A)') CLA8
                 CALL DZTEXT(MARE6Q,CHTEXT,0)
             ELSE
                 WRITE(CQINFO,'(I8)') LLINK
                 CALL DZTEXT(MARE7Q,CHTEXT,0)
             ENDIF
          ELSE
             IQUEST(10) = 0
          ENDIF
      ENDIF

  999 CQSTAK(MCQSIQ:) = CHSTAK
      RETURN
      END
