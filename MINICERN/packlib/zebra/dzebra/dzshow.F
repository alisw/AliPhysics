*
* $Id$
*
* $Log$
* Revision 1.2  1996/04/24 17:26:13  mclareni
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
      SUBROUTINE DZSHOW
     X (CHTEXT,IXSTOR,LBANK,CHOPT,ILINK1,ILINK2,IDATA1,IDATA2)
#include "zebra/bankparq.inc"
#include "zebra/questparq.inc"
#include "zebra/storparq.inc"
#include "zebra/mqsys.inc"
#include "zebra/qequ.inc"
#include "zebra/mzcn.inc"
#include "zebra/zbcdch.inc"
#include "zebra/zbcdk.inc"
#include "zebra/zunit.inc"
#include "zebra/dzc1.inc"
      CHARACTER *(*) CHOPT,CHTEXT

      PARAMETER ( NPDENQ =  3 )
      PARAMETER ( MPDCUQ =  1 )
      PARAMETER ( MPDNCQ =  2 )
      PARAMETER ( MPDNSQ =  3 )

      CHARACTER CHROUT*(*)
      PARAMETER (CHROUT = 'DZSHOW')

#include "zebra/q_jbit.inc"

#if (defined(CERNLIB_DEBUGON))&&(defined(CERNLIB_VFORT))
#include "zebra/debugvf2.inc"
#endif

      CQSTAK = CHROUT//'/'
      IQUEST(1) = 0

      IF (CHTEXT.NE.CDUMMQ)                                    THEN
          CQMAP(1)          = ' '
          CQMAP(2)(1:12)    = ' DZSHOW --- '
          CQMAP(2)(13:100)  = CHTEXT
          CQMAP(2)(101:110) = 'OPTIONS : '
          CQMAP(2)(111:130) = CHOPT
          CALL DZTEXT(0,CDUMMQ,2)
      ENDIF

      CALL DZOPT(CHOPT)

      CALL MZSDIV(IXSTOR,-1)


      IF (LBANK.EQ.0)                                      GO TO 999
      IF (JBIT(IQ(LBANK+KQS),IQDROP).EQ.1)                 GO TO 999


      CALL DZBKUP(0)
      IF (IQUEST(1).NE.0)                                  GO TO 999


      LWORK  = NQOFFS(1) + LQEND(1)  - NPDENQ - 1
      LWORKE = NQOFFS(1) + LQSTA(2)  - NPDENQ


      LCUR   = LBANK
      LEVEL  = 0
      LEVELO = -1


      MAXALL = 0
      NTBANK = 0
      DO 20 I=1,NDVMXQ
          IF(I.LE.JQDVLL.OR.I.GE.JQDVSY)
     I    MAXALL = MAXALL + LQEND(KQT+I) - LQSTA(KQT+I)
   20 CONTINUE



100   IF ( LEVEL.NE.LEVELO )  THEN


          IF(IFLOPT(MPOSVQ).NE.0) THEN
              CQMAP(1) = ' '
              WRITE(CQMAP(2),'('' DZSHOW  +++++ LEVEL '',I5,1X,
     W        10(''+''),T50,''Store  '',2A4,'' at absolute address '',
#if defined(CERNLIB_OCTMAP)
     W        O12,T100,10(''+'') )')
#endif
#if !defined(CERNLIB_OCTMAP)
     W        Z8,T100,10(''+'') )')
#endif
     W        LEVEL,NQSNAM(1),NQSNAM(2),
*  Map  addresses expressed in machine words
#if defined(CERNLIB_WORDMAP)
     W        LQSTOR+1
*  Map addresses expressed in bytes
#endif
#if !defined(CERNLIB_WORDMAP)
     W        (LQSTOR+1)*4
#endif
              CQMAP(3) = ' '
              CALL DZTEXT(0,CDUMMQ,3)
          ENDIF

          LEVELO = LEVEL

      ENDIF


      CALL DZSHPR (LCUR,ILINK1,ILINK2,IDATA1,IDATA2)
      IF (IQUEST(1).NE.0)                                  GO TO 999

      NTBANK = NTBANK + NL + ND + NBKOHQ
      IF (NTBANK.GE.MAXALL) THEN
          WRITE(CQINFO,'(I10,''/'',I10)') MAXALL,NTBANK
          CALL DZTEXT(MSHO1Q,CDUMMQ,0)
                                                           GO TO 999
      ENDIF


      LEVEL  = LEVEL+1
      LWORK  = LWORK + NPDENQ
      IF (LWORK.GE.LWORKE)             THEN
          WRITE(CQINFO,'(I10)') LEVEL
          CALL DZTEXT(MSHO2Q,CDUMMQ,0)
                                                           GO TO 999
      ENDIF


      IF (IFLOPT(MPOSVQ).EQ.0)   NS=0

      LQ(LWORK+MPDCUQ) = LCUR
      LQ(LWORK+MPDNCQ) = NS
      LQ(LWORK+MPDNSQ) = NS


200   IF (LQ(LWORK+MPDNCQ).LE.0)  THEN
          IF (LEVEL.GT.1.OR.IFLOPT(MPOSLQ).NE.0) THEN
              LCUR = LQ(KQS + LQ(LWORK+MPDCUQ))
              LEVEL  = LEVEL - 1
              LWORK  = LWORK - NPDENQ
              IF (LCUR.NE.LNULL)  THEN
                                                           GO TO 100
              ELSE
                  IF (LEVEL.GT.0) THEN
                                                           GO TO 200
                  ELSE
                                                           GO TO 999
                  ENDIF
              ENDIF
          ELSE
              LEVEL = LEVEL-1
              LWORK  = LWORK - NPDENQ
              IF (LEVEL.GT.0) THEN
                                                           GO TO 200
              ELSE
                                                           GO TO 999
              ENDIF
          ENDIF
      ENDIF


      LQ(LWORK+MPDNCQ) = LQ(LWORK+MPDNCQ) - 1
      LCUR = LQ(KQS+LQ(LWORK+MPDCUQ)
     X             -LQ(LWORK+MPDNSQ)+LQ(LWORK+MPDNCQ))

      IF (LCUR.EQ.LNULL)                                   GO TO 200
                                                           GO TO 100

  999 RETURN
      END
