*
* $Id$
*
* $Log$
* Revision 1.2  1996/04/24 17:26:06  mclareni
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
      SUBROUTINE DZCHST (CHTEXT,IXDIV,LBANK,CHOPT,ISUM)
      SAVE CHPART
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
      CHARACTER *(*) CHOPT,CHTEXT,CHPART(3)*6

      CHARACTER CHROUT*(*)
      PARAMETER (CHROUT = 'DZCHST')

      PARAMETER ( NPDENQ =  3 )
      PARAMETER ( MPDCUQ =  1 )
      PARAMETER ( MPDNCQ =  2 )
      PARAMETER ( MPDNSQ =  3 )

*       32 BIT MACHINES
#if defined(CERNLIB_B32)
      PARAMETER ( NBITS  =  8 )
*       36 BIT MACHINE
#endif
#if defined(CERNLIB_B36)
      PARAMETER ( NBITS  =  9 )
*       60 BIT MACHINE
#endif
#if defined(CERNLIB_B60)
      PARAMETER ( NBITS  = 15 )
*       64 BIT MACHINE
#endif
#if defined(CERNLIB_B64)
      PARAMETER ( NBITS  = 16 )
#endif

      PARAMETER ( NLSUMQ = 6 )
      INTEGER ISUM(*),ISMOLD(NLSUMQ)

      DATA CHPART /'DATA','LINK','SYSTEM'/

#if (defined(CERNLIB_DEBUGON))&&(defined(CERNLIB_VFORT))
#include "zebra/debugvf2.inc"
#endif

      CQSTAK = CHROUT//'/'
      IQUEST(1) = 0

      CALL DZOPT(CHOPT)

      IF (IFLOPT(MPOSVQ).NE.0)  THEN
          DO 10 I=1,NLSUMQ
  10      ISMOLD(I) = ISUM(I)
      ENDIF

      CALL MZSDIV(IXDIV,-1)


      LWORK  = NQOFFS(1) + LQEND(1)  - NPDENQ - 1
      LWORKE = NQOFFS(1) + LQSTA(2)  - NPDENQ


      LCUR   = LBANK
      LEVEL  = 0


      MAXALL = 0
      NTBANK = 0
      DO 20 I=1,NDVMXQ
          IF(I.LE.JQDVLL.OR.I.GE.JQDVSY)
     I    MAXALL = MAXALL + LQEND(KQT+I) - LQSTA(KQT+I)
   20 CONTINUE



      CALL VZERO(IQUEST(71),12)


  100 CALL MZCHLS (NCHEKQ,LCUR)
      IF (IQFOUL.NE.0) THEN
          CALL DZBKDV(LCUR)
          IF (IQUEST(1).NE.0)                              GO TO 999
          WRITE(CQINFO,'(A,''/'',I8)') CQDIV,LCUR
          CALL DZTEXT(MSHP1Q,CDUMMQ,0)
          IQUEST(1) = 1
                                                           GO TO 999
      ENDIF

      CALL DZCHV1 (LCUR+1+NOFLIQ+KQS,LCUR+IQND+NOFLIQ+KQS,70,0)
      IF (IQUEST(1).NE.0)                                  GO TO 999
      CALL DZCHV1 (LCUR-IQNL+KQS,LCUR+2+KQS,74,0)
      IF (IQUEST(1).NE.0)                                  GO TO 999
      CALL DZCHV1 (LCUR+3+KQS,LCUR+NOFLIQ+KQS,78,0)
      IF (IQUEST(1).NE.0)                                  GO TO 999
      CALL DZCHV1 (LCUR-IQNL-1+KQS,LCUR-IQNL-1+KQS,78,0)
      IF (IQUEST(1).NE.0)                                  GO TO 999

      NTBANK = NTBANK + NL + ND + NBKOHQ
      IF (NTBANK.GE.MAXALL) THEN
          WRITE(CQINFO,'(I10,''/'',I10)') MAXALL,NTBANK
          CALL DZTEXT(MSHO1Q,CDUMMQ,0)
          IQUEST(1) = 1
                                                           GO TO 999
      ENDIF


      LEVEL  = LEVEL+1
      LWORK  = LWORK + NPDENQ
      IF (LWORK.GE.LWORKE)             THEN
          WRITE(CQINFO,'(I10)') LEVEL
          CALL DZTEXT(MSHO2Q,CDUMMQ,0)
          IQUEST(1) = 1
                                                           GO TO 999
      ENDIF


      IF (IFLOPT(MPOSDQ).EQ.0)   IQNS = 0

      LQ(LWORK+MPDCUQ) = LCUR
      LQ(LWORK+MPDNCQ) = IQNS
      LQ(LWORK+MPDNSQ) = IQNS


  200 IF (LQ(LWORK+MPDNCQ).LE.0)  THEN
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
                                                           GO TO 300
                  ENDIF
              ENDIF
          ELSE
              LEVEL = LEVEL-1
              LWORK  = LWORK - NPDENQ
              IF (LEVEL.GT.0) THEN
                                                           GO TO 200
              ELSE
                                                           GO TO 300
              ENDIF
          ENDIF
      ENDIF


      LQ(LWORK+MPDNCQ) = LQ(LWORK+MPDNCQ) - 1
      LCUR = LQ(KQS+LQ(LWORK+MPDCUQ)
     X             -LQ(LWORK+MPDNSQ)+LQ(LWORK+MPDNCQ))

      IF (LCUR.EQ.LNULL)                                   GO TO 200
                                                           GO TO 100



  300 DO 320 I=1,NLSUMQ/2
          II = (I-1)*2

          JCARRY = 0
          IBIT   = 1
          DO 310 JF=70+I*4-3,70+I*4
              IQUEST (JF) = IQUEST(JF) + JCARRY
              JCARRY = IQUEST(JF)/2**NBITS
              IQUEST (JF) = IQUEST(JF) - JCARRY*(2**NBITS)
              CALL SBYT(IQUEST(JF),ISUM(II+1),IBIT,NBITS)
  310     IBIT   = IBIT + NBITS
          ISUM(II+2) = JCARRY
  320 CONTINUE

      IF (IFLOPT(MPOSVQ).NE.0) THEN
          DO 400 I=1,NLSUMQ/2
              II = (I-1)*2 + 1
              IF (     (ISUM(II)  .NE.ISMOLD(II)  )
     I            .OR. (ISUM(II+1).NE.ISMOLD(II+1)) ) THEN
                  IQUEST(10+I) = 1
                  IQUEST(1)    = 1
              ELSE
                  IQUEST(10+I) = 0
              ENDIF
  400     CONTINUE
      ENDIF

      IF (CHTEXT.NE.CDUMMQ)     THEN
          CQMAP(1)        = ' '
          CQMAP(2)(1:10)  = ' * '//CHROUT//' '
          CQMAP(2)(11:29) = CHTEXT
          CQMAP(2)(30:41) = ' / OPTION : '
          CQMAP(2)(42:47) = CHOPT
          DO 500 I=1,NLSUMQ/2
              II = (I-1)*2 + 1
              IF (IFLOPT(MPOSVQ).NE.0) THEN
                  IF (IQUEST(10+I).NE.0)   THEN
                      CQMAP(2)(118:130) = '??PROBLEMS? '
                  ELSE
                      CQMAP(2)(118:130) = '   OK'
                  ENDIF
                  WRITE(CQMAP(2)(48:117),
     W             '(''OLD='',Z4,1X,Z16,4X,''NEW='',Z4,1X,Z16,T63,A)')
     W             ISMOLD(II+1),ISMOLD(II),ISUM(II+1),ISUM(II),CHPART(I)
              ELSE

                  WRITE(CQMAP(2)(48:),'(''NEW='',Z4,1X,Z16,T63,A)')
     W             ISUM(II+1),ISUM(II),CHPART(I)
              ENDIF
              CALL DZTEXT(0,CDUMMQ,2)
              CQMAP(2) = ' '
  500     CONTINUE
      ENDIF

  999 RETURN
      END
