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
*-----------------------------------------------------------
#include "zebra/pilot.h"
#if (defined(CERNLIB_DEBUGON))&&(defined(CERNLIB_VFORT))
#include "zebra/debugvf1.inc"
#endif
      SUBROUTINE DZSURV (CHTEXT,IXDIV,LBANK)
      SAVE NWDES
#include "zebra/bankparq.inc"
#include "zebra/divparq.inc"
#include "zebra/questparq.inc"
#include "zebra/mqsys.inc"
#include "zebra/qequ.inc"
#include "zebra/mzcn.inc"
#include "zebra/zbcd.inc"
#include "zebra/zbcdk.inc"
#include "zebra/zunit.inc"
#include "zebra/dzc1.inc"
      DIMENSION       NEWID(10)

      CHARACTER CHTEXT*(*)

      CHARACTER CHROUT*(*)
      PARAMETER (CHROUT = 'DZSURV')

#include "zebra/q_jbit.inc"

#if (defined(CERNLIB_DEBUGON))&&(defined(CERNLIB_VFORT))
#include "zebra/debugvf2.inc"
#endif

      DATA  NWDES  / 4 /


      CQSTAK = CHROUT//'/'
      IQUEST(1) = 0

      CALL MZSDIV (IXDIV,-1)


      IF (LBANK.EQ.0)                                      GO TO 999
      IF (JBIT(IQ(LBANK+KQS),IQDROP).EQ.1)                 GO TO 999

      LSTART = LBANK

      CALL MZCHLS(NCHEKQ,LSTART)
      IF (IQFOUL.NE.0)                                     GO TO 91
      IDSTR  = IQID


      IF (CHTEXT.NE.CDUMMQ.AND.CHTEXT.NE.'-')        THEN
          CALL ZPAGE(IQPRNT,7)
          CQMAP(1)        = ' '
          CQMAP(2)        = ' DZSURV --- '
          CQMAP(2)(13:99) = CHTEXT
          WRITE(CQMAP(2)(100:),'('' ST= '',2A4,''  LSTART= '',I8)')
     X      NQPNAM(KQT+1),NQPNAM(KQT+2),LSTART
          CALL DZTEXT(0,CDUMMQ,2)
      ENDIF

      IF (CHTEXT.NE.'-') THEN
          CALL ZPAGE(IQPRNT,5)
          CQMAP(1) = ' '
          CQMAP(2) = '   NWCUM     NW   WBK  NBK    IDENTIFIER(S)'
          CQMAP(3) = ' '
          CALL DZTEXT(0,CDUMMQ,3)
      ENDIF

      MAXALL = 0
      DO 20 I=1,NDVMXQ
          IF(I.LE.JQDVLL.OR.I.GE.JQDVSY)
     I    MAXALL = MAXALL + LQEND(KQT+I) - LQSTA(KQT+I)
   20 CONTINUE
      NBCUM    = 0
      NWCUM    = 0
      NEWNW    = 0
      NEWND    = 0
      NEWNBK   = 0
      NEWRIN   = IQBLAN
      NEWLOW   = 99
      NEWUP    = 1
      NEWNID   = 0
      NEWID(1) = 0
      NEWLEV   = 1
      LEVELH   = 0
      LTAB1    = NQOFFT(1)
      LEND1    = LQEND(LTAB1+1)
      LWORK    = LQWKTB
      NWORK    = NQWKTB
      NLVMAX   = NWORK/NWDES
      LDESC    = LWORK - NWDES
      LKBCD    = IQBLAN
      LGO      = LSTART
      L        = LSTART

                                                           GO TO 52
   31 NEWLEV = NEWLEV - 1

      IF (CHTEXT.NE.'-') CALL DZTEXT(1,CDUMMQ,1)

      IF (NEWLEV.EQ.1)                                     GO TO 9991
      LDESC = LDESC - NWDES
   34 IF (LQ(LDESC).EQ.LQ(LDESC+1))                        GO TO 31
      LQ(LDESC) = LQ(LDESC) + 1
   39 CALL VBLANK (IQUEST(2),3)
      CALL USET   (LQ(LDESC),IQUEST,2,4)
      IQUEST(1) = IQMINS
      CALL URIGHT (IQUEST,1,4)
      CALL UBUNCH (IQUEST,LKBCD,4)
      NEWNW    = 0
      NEWND    = 0
      NEWNBK   = 0
      NEWRIN   = IQBLAN
      NEWLOW   = 99
      NEWUP    = 1
      NEWNID   = 0
      NEWID(1) = 0
      LEVELH   = 1
      LDESH    = LWORK
      LGO      = LSTART
                                                           GO TO 42
   41 LEVELH      = LEVELH + 1
      LDESH       = LDESH  + NWDES
   42 JFOLL       = LQ(LDESH)
      LQ(LDESH+2) = LGO
      LQ(LDESH+3) = LGO
      LM          = LGO
                                                           GO TO 46
   43 NEWRIN = IQLETT(MPOSRQ)
   44 IF (LEVELH.EQ.0)                                     GO TO 71
   45 LM = LQ(LM+KQS)
      IF (LM.EQ.0)                                         GO TO 47
      IF (LM.EQ.LQ(LDESH+3))                               GO TO 47
   46 CALL MZCHLS (NCHEKQ,LM)
      IF (IQFOUL.NE.0)                                     GO TO 91
      IF (IQNS-JFOLL.LT.0)                                 GO TO 45
      LQ(LDESH+2) = LM
      K   = LM - JFOLL
      LGO = LQ(K+KQS)

      IF (LGO.EQ.0)                                        GO TO 45
      IF (LEVELH+1.LT.NEWLEV)                              GO TO 41
      L = LGO
                                                           GO TO 52
   47 LEVELH = LEVELH - 1
      IF (LEVELH.EQ.0)                                     GO TO 71
      LDESH = LDESH - NWDES
      JFOLL = LQ(LDESH)
      LM    = LQ(LDESH+2)
                                                           GO TO 45
   51 K = L
      L = LQ(K+KQS)
      IF (L.EQ.0)                                          GO TO 44
      IF (L.EQ.LGO)                                        GO TO 43
   52 CALL MZCHLS(NCHEKQ,L)
      IF (IQFOUL.NE.0)                                     GO TO 91
      NEWNBK = NEWNBK + 1
      N      = NBKOHQ + IQNL + IQND
      NEWND  = MAX(NEWND,N)
      NEWNW  = NEWNW + N
      IF (NEWNW.GE.MAXALL)                                 GO TO 71
      IF (IQID.EQ.NEWID(1))                                GO TO 57
      IF (NEWNID.EQ.0)                                     GO TO 56
      IF (NEWNID.EQ.10)                                    GO TO 57
      IF (IUCOMP(IQID,NEWID,NEWNID).NE.0)                  GO TO 57
   56 NEWNID = NEWNID + 1
      NEWID(NEWNID) = IQID
   57 IF (IQNS.EQ.0)                                       GO TO 51
      N = MIN(NEWLOW-1,IQNS)
      DO 62 J=1,N
          IF (LQ(L-J+KQS).NE.0)        THEN
              NEWLOW = J
                                                           GO TO 64
          ENDIF
   62 CONTINUE
   64 JA = MAX(NEWLOW,NEWUP) + 1
      DO 65 J=JA,IQNS
          IF (LQ(L-J+KQS).NE.0)  NEWUP=J
   65 CONTINUE

                                                           GO TO 51
   71 IF (NEWNBK.EQ.0)                                     GO TO 75
      NBCUM = NBCUM + NEWNBK
      NWCUM = NWCUM + NEWNW
      CALL VBLANK(IQUEST,NEWLEV)
      IQUEST(NEWLEV) = LKBCD
      NEWNID = MIN(NEWNID,21-NEWLEV)
      IF (CHTEXT.NE.'-') THEN
          WRITE(CQLINE,'(1X,2I7,I6,I5,2(1X,A1),20(1X,A4))')
     W                   NWCUM,NEWNW,NEWND,NEWNBK,NEWRIN
     W,                  (IQUEST(J),J=1,NEWLEV)
     W,                  (NEWID(J), J=1,NEWNID)
          CALL DZTEXT(0,CDUMMQ,1)
      ENDIF

      IF (NWCUM.GT.MAXALL) THEN
          WRITE(CQINFO,'(I10,''/'',I10)') NWCUM,MAXALL
          CALL DZTEXT(MSHO1Q,CDUMMQ,0)
          IQUEST(1) = 1
                                                           GO TO 999
      ENDIF

      IF (NEWLOW.GE.64)                                    GO TO 75
      IF (NEWLEV.GE.NLVMAX)                                GO TO 74
      NEWLEV = NEWLEV + 1
      LDESC  = LDESC  + NWDES
      LQ(LDESC)   = NEWLOW
      LQ(LDESC+1) = MAX(NEWLOW,NEWUP)
                                                           GO TO 39
   74 CALL DZTEXT(MSUR1Q,CDUMMQ,0)
   75 IF (NEWLEV.NE.1)                                     GO TO 34
                                                           GO TO 9991
   91 IQUEST(1) = LGO
      IQUEST(2) = K
      IQUEST(3) = L
      CALL ZFATAM (CHROUT)


 9991 IF (CHTEXT.NE.'-') THEN
          WRITE(CQMAP,'(1X,/,'' DZSURV --- The structure supported by'',
     +    '' bank '',A4,'' at '',I10,'' in store '',2A4,'' occupies '',
     +    I10,'' words in '',I6,'' banks '')')
     +    IDSTR,LSTART,NQPNAM(KQT+1),NQPNAM(KQT+2),NWCUM,NBCUM
          CALL DZTEXT(0,CDUMMQ,2)
      ENDIF
      IQUEST(11) = IDSTR
      IQUEST(12) = LSTART
      IQUEST(13) = NQPNAM(KQT+1)
      IQUEST(14) = NQPNAM(KQT+2)
      IQUEST(15) = NWCUM
      IQUEST(16) = NBCUM
  999 RETURN
      END
