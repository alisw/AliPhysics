*
* $Id$
*
* $Log$
* Revision 1.2  1996/04/18 16:11:58  mclareni
* Incorporate changes from J.Zoll for version 3.77
*
* Revision 1.1.1.1  1996/03/06 10:47:19  mclareni
* Zebra
*
*
#include "zebra/pilot.h"
      SUBROUTINE MZPUSH (IXDIV,LORGP,INCNLP,INCNDP,CHOPT)

C-    Change the size of a bank, user called

#include "zebra/zbcd.inc"
#include "zebra/zstate.inc"
#include "zebra/zunit.inc"
#include "zebra/zvfaut.inc"
#include "zebra/mqsys.inc"
#include "zebra/mzcl.inc"
#include "zebra/mzcn.inc"
#include "zebra/mzct.inc"
C--------------    End CDE                             --------------
      DIMENSION    IXDIV(9),LORGP(9),INCNLP(9),INCNDP(9)
      CHARACTER    *(*) CHOPT
#if (defined(CERNLIB_QTRHOLL))&&(!defined(CERNLIB_A6M))
      DIMENSION    NAMESR(2)
      DATA  NAMESR / 4HMZPU, 4HSH   /
#endif
#if (defined(CERNLIB_QTRHOLL))&&(defined(CERNLIB_A6M))
      DATA  NAMESR / 6HMZPUSH /
#endif
#if !defined(CERNLIB_QTRHOLL)
      CHARACTER    NAMESR*8
      PARAMETER   (NAMESR = 'MZPUSH  ')
#endif

#include "zebra/q_jbit.inc"
#include "zebra/q_jbyt.inc"
#include "zebra/q_sbit1.inc"
#include "zebra/q_sbyt.inc"

#include "zebra/qtrace.inc"

      IF (IXDIV(1).EQ.-7)          GO TO 12
      CALL MZSDIV (IXDIV,0)

   12 CALL MZCHNB (LORGP)
      LORG  = LORGP(1)
      INCNL = INCNLP(1)
      INCND = INCNDP(1)

      CALL UOPTC (CHOPT,'RI',IQUEST)
      IFLAG = MIN (2, IQUEST(1)+2*IQUEST(2))

C-         IFLAG = 0  general
C-                 1  R-educe
C-                 2  I-solated

#if defined(CERNLIB_QDEVZE)
      IF (NQDEVZ.NE.0)
     +WRITE (IQLOG,9809) LORG,INCNL,INCND,IFLAG
 9809 FORMAT (1X/' DEVZE MZPUSH,  Entry for LORG,INCNL,INCND,IFLAG= '
     F,5I8)
#endif
#if defined(CERNLIB_QDEBUG)
      IF (IQVSTA.NE.0)       CALL ZVAUTX
#endif
      IF ((INCNL.EQ.0) .AND. (INCND.EQ.0))  GO TO 999
      LQSYSR(KQT+1) = LORG

C--                Find division

      JQDIVI = MZFDIV (-7, LORG)
      IF (JQDIVI.EQ.0)             GO TO 91

C--                Set bank parameters

#if defined(CERNLIB_QDEBUG)
      CALL MZCHLS (-7,LORG)
      IF (IQFOUL.NE.0)             GO TO 91
      NL    = IQNL
      NS    = IQNS
      ND    = IQND
      NQNIO = IQNIO
      NQID  = IQID
#endif
#if defined(CERNLIB_QDEVZE)
      IF (NQDEVZ.GE.7)
     +WRITE (IQLOG,9831) JQDIVI,IQLS,IQID,IQNL,IQNS,IQND
 9831 FORMAT (16X,'JQDIVI,IQLS,IQID,IQNL,IQNS,IQND=',I3,I7,1X,A4,4I8)
#endif
#if !defined(CERNLIB_QDEBUG)
      NQID  = IQ(KQS+LORG-4)
      NL    = IQ(KQS+LORG-3)
      NS    = IQ(KQS+LORG-2)
      ND    = IQ(KQS+LORG-1)
      NQNIO = JBYT (IQ(KQS+LORG),19,4)
#endif
      NQNL = NL + INCNL
      NQNS = MIN (NS,NQNL)
      NQND = ND + INCND
      IF (NS.EQ.NL)  NQNS = NQNL
#if defined(CERNLIB_QDEBPRI)
      IF (NQLOGL.GE.2)
     + WRITE (IQLOG,9032) JQSTOR,JQDIVI,LORG,NQID,INCNL,INCND,CHOPT
 9032 FORMAT (' MZPUSH-  Store/Div',2I3,' L/ID/INCNL/INCND/OPT=',
     FI9,1X,A4,2I7,1X,A)
#endif
      IF (JBIT(IQ(KQS+LORG),IQDROP).NE.0)   GO TO 92

C--                Check for bad parameters

      IF (NQND+NQNL.GE.LQSTA(KQT+21)) GO TO 93
      IF (NQND.LT.0)               GO TO 93
      IF (NQNL.GT.64000)           GO TO 93
      IF (NQNS.LT.0)               GO TO 93

      NLC = MIN (NL,NQNL)
      NSC = MIN (NS,NQNS)
      NDC = MIN (ND,NQND)

C--                Check giving up non-zero structural links

      IF (NQNS.GE.NS)              GO TO 36
      L  = LORG - NS - 1
      LD = LORG - NQNS
   34 L  = L + 1
      IF (L.GE.LD)                 GO TO 36
      LNZ = LQ(KQS+L)
   35 IF (LNZ.EQ.0)                GO TO 34
      IF (LQ(KQS+LNZ+2).NE.L)         GO TO 34
      IF (JBIT(IQ(KQS+LNZ),IQDROP).EQ.0)   GO TO 94
      LNZ = LQ(KQS+LNZ)
      GO TO 35

C--                Ready I/O characteristic

   36 LN = LORG - NL - NQNIO - 1
      CALL UCOPY (LQ(KQS+LN),NQIOCH,NQNIO+1)
      IF (NQNIO.NE.0)  NQIOSV(1)=0
      NQIOCH(1) = MSBYT (NQNL+NQNIO+12,NQIOCH(1),1,16)

C--                Re-enter after garbage collection, if any

   41 LE = LORG + ND + 9

C------            Check for short-cuts

      INCTT = INCNL + INCND
      INCMX = MAX (INCNL,INCND)
      INCMI = MIN (INCNL,INCND)
      CALL MZRESV
      IF (JQMODE.NE.0)             GO TO 45

C--                Last bank in forward division

      IF (LE.NE.LQEND(KQT+JQDIVI))    GO TO 51
      IF (INCNL.GE.0)              GO TO 42
      IF (IFLAG.NE.1)              GO TO 42
      IF ((NQRESV.GE.INCTT).AND.(NQRESV.LT.INCND))  GO TO 42

      LNN = LN - INCNL
      CALL UCOPY (NQIOCH,LQ(KQS+LNN),NQNIO+1)
      IQ(KQS+LORG-3) = NQNL
      IQ(KQS+LORG-2) = NQNS

      NWD = -INCNL
      CALL MZPUDX (LN,NWD)
      INCNL = 0
      INCTT = INCND
      LN    = LNN
      NL    = NQNL

   42 NQRESV = NQRESV - INCTT
      IF (NQRESV.LT.0)             GO TO 49
      NDELTA = INCNL
      LNEW   = LORG + NDELTA
      LQEND(KQT+JQDIVI) = LQEND(KQT+JQDIVI) + INCTT
#if defined(CERNLIB_QDEVZE)
      IF (NQDEVZ.GE.7)  WRITE (IQLOG,9848) NDELTA,LNEW
#endif
      IF (NDELTA.EQ.0)  THEN
          IQ(KQS+LNEW-1) = NQND
          IF (IFLAG.NE.0)          GO TO 81
          IF (INCMI.GE.0)          GO TO 81
          GO TO 71
        ELSE
          CALL UCOPY2 (LQ(KQS+LORG-NLC),LQ(KQS+LNEW-NLC),NLC+NDC+9)
          IF (INCNL.GT.0)  CALL VZERO (LQ(KQS+LNEW-NQNL),INCNL)
          LQ(KQS+LN)     = NQIOCH(1)
          IQ(KQS+LNEW-3) = NQNL
          IQ(KQS+LNEW-2) = NQNS
          IQ(KQS+LNEW-1) = NQND
          GO TO 61
        ENDIF

C--                First bank in reverse division

   45 IF (LN.NE.LQSTA(KQT+JQDIVI))    GO TO 51
      IF (INCND.GE.0)              GO TO 47
      IF (IFLAG.NE.1)              GO TO 47
      IF ((NQRESV.GE.INCTT).AND.(NQRESV.LT.INCNL))  GO TO 47
      IQ(KQS+LORG-1) = NQND

      L   = LE + INCND
      NWD = -INCND
      CALL MZPUDX (L,NWD)
      INCND = 0
      INCTT = INCNL
      ND    = NQND

   47 NQRESV = NQRESV - INCTT
      IF (NQRESV.LT.0)             GO TO 49
      LNN    = LN - INCTT
      NDELTA = -INCND
      LQSTA(KQT+JQDIVI) = LNN

      LNEW  = LORG + NDELTA

      IF (NDELTA.NE.0)  CALL UCOPY2 (LQ(KQS+LORG-NLC)
     +,                              LQ(KQS+LNEW-NLC), NLC+NDC+9)

      IF (INCNL.GT.0)  CALL VZERO (LQ(KQS+LNEW-NQNL),INCNL)
      CALL UCOPY (NQIOCH,LQ(KQS+LNN),NQNIO+1)

      IQ(KQS+LNEW-3) = NQNL
      IQ(KQS+LNEW-2) = NQNS
      IQ(KQS+LNEW-1) = NQND
#if defined(CERNLIB_QDEVZE)
      IF (NQDEVZ.GE.7)  WRITE (IQLOG,9848) NDELTA,LNEW
 9848 FORMAT (' DEVZE MZPUSH,  Edge bank with NDELTA,LNEW=',2I8)
#endif

      IF (NDELTA.NE.0)             GO TO 61
      IF (IFLAG.NE.0)              GO TO 81
      IF (INCMI.GE.0)              GO TO 81
      GO TO 71

C--                Garbage collection

   49 CALL MZGAR1
      LORG = LQSYSR(KQT+1)
      LN   = LORG - NL - NQNIO - 1
#if defined(CERNLIB_QDEVZE)
      IF (NQDEVZ.GE.7)
     +WRITE (IQLOG,9849)
 9849 FORMAT (1X/' DEVZE MZPUSH,  Garbage collected for edge bank')
#endif
      GO TO 41

C----              Reduction only

   51 IF (INCMX.GT.0)              GO TO 56
      IF (INCNL.EQ.0)              GO TO 52

C--                Link part

      LNN = LN - INCNL
      CALL UCOPY (NQIOCH,LQ(KQS+LNN),NQNIO+1)
      IQ(KQS+LORG-3)= NQNL
      IQ(KQS+LORG-2)= NQNS

      CALL MZPUDX (LN,-INCNL)

#if defined(CERNLIB_QDEVZE)
      IF (NQDEVZ.GE.7)
     +WRITE (IQLOG,9851) INCNL
 9851 FORMAT (' DEVZE MZPUSH,  In-situ links with INCNL=',I8)
#endif
      IF (INCND.EQ.0)              GO TO 54

C--                Data part

   52 IQ(KQS+LORG-1) = NQND
      LD  = LE + INCND
      NWD = -INCND
      CALL MZPUDX (LD,NWD)
#if defined(CERNLIB_QDEVZE)
      IF (NQDEVZ.GE.7)
     +WRITE (IQLOG,9852) INCND
 9852 FORMAT (' DEVZE MZPUSH,  In-situ data with INCND=',I8)
#endif

   54 LNEW   = LORG
      NDELTA = 0
      IF (IFLAG.NE.0)              GO TO 999
      GO TO 71

C------            Lift replacement bank

   56 J = 64*(32*NQNIO + NQNIO + 1) + 1
      NQIOCH(1) = MSBYT (J,NQIOCH(1),1,16)

      NQBIA = 2
      CALL MZLIFT (-7,LNEW,0,63,NQID,-1)
      LORG   = LQSYSR(KQT+1)
      NDELTA = LNEW - LORG

      CALL UCOPY (LQ(KQS+LORG-NLC),LQ(KQS+LNEW-NLC),NLC+4)
      CALL UCOPY (IQ(KQS+LORG),    IQ(KQS+LNEW),    NDC+1)
      IQ(KQS+LORG) = MSBIT1 (IQ(KQS+LORG),IQDROP)
#if defined(CERNLIB_QDEVZE)
      IF (NQDEVZ.GE.7)
     +WRITE (IQLOG,9857) LORG,LNEW
 9857 FORMAT (' DEVZE MZPUSH,  Push by copy LORG -> LNEW=',2I8)
#endif

C------            Up-date immediate links only

   61 IF (IFLAG.LT.2)              GO TO 71
#if defined(CERNLIB_QDEVZE)
      IF (NQDEVZ.GE.7)
     +WRITE (IQLOG,9861)
 9861 FORMAT (' DEVZE MZPUSH,  Update immediate links only')
#endif

C----              Update according to k-link in pushed bank

      K = LQ(KQS+LNEW+2)
      IF (K.EQ.0)                  GO TO 62
      IF (LQ(KQS+K).NE.LORG)          GO TO 95
      LQ(KQS+K) = LNEW

C----              Update according to link 0

   62 K = LNEW
      L = LQ(KQS+K)
      IF (L.EQ.0)                  GO TO 65
      IF (L.EQ.LORG)               GO TO 64
      LQ(KQS+L+2) = K
   63 K = L
      L = LQ(KQS+K)
      IF (L.EQ.0)                  GO TO 65
      IF (L.NE.LORG)               GO TO 63
   64 LQ(KQS+K) = LNEW

C----              Update k- and up-link in vertically dependent banks

   65 K = LNEW - NSC - 1

C--                          each link
   66 K = K + 1
      IF (K.GE.LNEW)               GO TO 81
      L = LQ(KQS+K)
      IF (L.EQ.0)                  GO TO 66
      IF (LQ(KQS+L+2).NE.K-NDELTA)    GO TO 66
      LQ(KQS+L+2) = K

C--                          and its linear structure

      LF = L
   68 LQ(KQS+L+1) = LNEW
      L = LQ(KQS+L)
      IF (L.EQ.LF)                 GO TO 66
      IF (L.NE.0)                  GO TO 68
      GO TO 66

C------            Global update of links

   71 MQDVGA = 0
      MQDVWI = 0
      JQSTMV = -1
#if defined(CERNLIB_QDEBPRI)
      IF (NQLOGL.GE.1)
     + WRITE (IQLOG,9071) JQSTOR,JQDIVI,LORG,NQID
 9071 FORMAT (' MZPUSH-  Store/Div',2I3,' Relocation pass for L/ID ='
     F,I9,1X,A4)
#endif
#if defined(CERNLIB_QDEVZE)
      IF (NQDEVZ.GE.7)
     +WRITE (IQLOG,9871)
 9871 FORMAT (' DEVZE MZPUSH,  Update by relocation pass')
#endif

      CALL MZTABM

      LMT  = LQMTA - 8
   74 LMT  = LMT + 8
      IF (LQ(LMT).NE.JQDIVI)       GO TO 74
      LQ(LMT+1) = 2

      CALL MZTABX
      LQMTE = LQMTLU

      LQ(LQTA-1) = LORG - NL - NQNIO - 1
      LQ(LQTA)   = LORG - NLC
      LQ(LQTA+1) = LORG + NDC + 9
      LQ(LQTA+2) = NDELTA
      LQ(LQTA+3) = 0
      LQ(LQTA+4) = LORG + ND + 9

      LQTE  = LQTA + 4

      CALL MZRELX

      NQDPSH(KQT+JQDIVI) = NQDPSH(KQT+JQDIVI) + 1

C------            Finished, reset LORG, clear new data words

   81 LORGP(1) = LNEW
      IF (INCND.GT.0)  CALL VZERO (IQ(KQS+LNEW+ND+1),INCND)
#include "zebra/qtrace99.inc"
      RETURN

C----              Error conditions

   95 NQCASE = 3
      NQFATA = 1
      IQUEST(19) = K
      GO TO 92

   94 NQCASE = 1
      NQFATA = 2
      IQUEST(19) = L - LORG
      IQUEST(20) = LQ(KQS+L)
   93 NQCASE = NQCASE + 1
   92 NQCASE = NQCASE + 1
      NQFATA = NQFATA + 7
      IQUEST(12) = NQID
      IQUEST(13) = NS
      IQUEST(14) = NL
      IQUEST(15) = ND
      IQUEST(16) = NQNIO
      IQUEST(17) = INCNL
      IQUEST(18) = INCND
   91 NQCASE = NQCASE + 1
      NQFATA = NQFATA + 1
      IQUEST(11) = LORG
#include "zebra/qtofatal.inc"
      END
*      ==================================================
#include "zebra/qcardl.inc"
