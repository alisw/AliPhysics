*
* $Id$
*
* $Log$
* Revision 1.2  1996/04/18 16:11:37  mclareni
* Incorporate changes from J.Zoll for version 3.77
*
* Revision 1.1.1.1  1996/03/06 10:47:18  mclareni
* Zebra
*
*
#include "zebra/pilot.h"
      SUBROUTINE MZLIFT (IXDIV,LP,LSUPP,JBIAS,NAME,NZERO)

C-    Lift a bank, user called

#include "zebra/zmach.inc"
#include "zebra/zstate.inc"
#include "zebra/zunit.inc"
#include "zebra/zvfaut.inc"
#include "zebra/mqsys.inc"
#include "zebra/eqlqmst.inc"
#include "zebra/eqlqform.inc"
#include "zebra/mzcl.inc"
#include "zebra/mzcn.inc"
#include "zebra/mzct.inc"
C--------------    End CDE                             --------------
      DIMENSION    IXDIV(9), LP(9), LSUPP(9), NAME(9)
#if (defined(CERNLIB_QTRHOLL))&&(!defined(CERNLIB_A6M))
      DIMENSION    NAMESR(2)
      DATA  NAMESR / 4HMZLI, 4HFT   /
#endif
#if (defined(CERNLIB_QTRHOLL))&&(defined(CERNLIB_A6M))
      DATA  NAMESR / 6HMZLIFT /
#endif
#if !defined(CERNLIB_QTRHOLL)
      CHARACTER    NAMESR*8
      PARAMETER   (NAMESR = 'MZLIFT  ')
#endif

#include "zebra/q_jbit.inc"
#include "zebra/q_jbyt.inc"
#include "zebra/q_sbyt.inc"
#include "zebra/q_shiftl.inc"
#include "zebra/q_locf.inc"


#include "zebra/qtrace.inc"
#if defined(CERNLIB_QDEBUG)
      IF (IQVSTA.NE.0)       CALL ZVAUTX
#endif

C--                Copy parameters

      IF (JBIAS.NE.63)  THEN
          NQBIA = JBIAS
          NIO   = JBYT (NAME(5),12,4)
          CALL UCOPY (NAME,NQID,NIO+5)
          IF (NIO.NE.0)  NQIOSV(1)=0
        ENDIF

      JDV   = IXDIV(1)
      LQSUP = LSUPP(1)
      IF (NQBIA.GE.2)  LQSUP = 0
      ICHORG = NQIOCH(1)
      NTOT   = NQNL + NQND + 10

      IF (JDV.EQ.-7)                 GO TO 24
      IF (JBYT(JDV,27,6).NE.JQSTOR)  GO TO 22
      JQDIVI = JBYT (JDV,1,26)
      IF (JQDIVI.LT.21)              GO TO 23
   22 CALL MZSDIV (JDV,0)
   23 CALL MZCHNB (LP)
   24 CONTINUE

C--                Check  bank parameters

#if defined(CERNLIB_QDEBUG)
      J = JBYT (NQID,IQBITW-7,8)
      IF (J.EQ.0)                  GO TO 91
      IF (J.EQ.255)                GO TO 91
      IF (NTOT.GE.LQSTA(KQT+21))      GO TO 91
      IF (NQNS.GT.NQNL)            GO TO 91
      IF (NQNS.LT.0)               GO TO 91
      IF (NQNL.GT.64000)           GO TO 91
      IF (NQND.LT.0)               GO TO 91
      IF (NQBIA.GE.3)              GO TO 91

      IF (LQSUP.EQ.0)              GO TO 25
      CALL MZCHLS (-7,LQSUP)
      IF (IQFOUL.NE.0)             GO TO 92
      IF (NQBIA.EQ.1)              GO TO 26
      IF (JBIT(IQ(KQS+LQSUP),IQDROP).NE.0)  GO TO 92
      IF (IQNS+NQBIA.LT.0)         GO TO 93
      GO TO 26

   25 IF (NQBIA.LE.0)              GO TO 92
   26 CONTINUE
#endif

C----       Find   LNEXT,  future 'next' bank
C-                 LSAME,  a bank in the same linear structure
C-                    LS,  division selecting bank
C-                   IDN,  numeric ID

      IDN   = 1
      LS    = LQSUP
      LSAME = LQSUP
      LNEXT = LQSUP
      IF (NQBIA.GT.0)              GO TO 38
      LNEXT = LQ(KQS+LNEXT+NQBIA)
      IF (LNEXT.EQ.0)              GO TO 36
      LSAME = LNEXT
      LS    = LNEXT
#if defined(CERNLIB_QDEBUG)
      CALL MZCHLS (-7,LNEXT)
      IF (IQFOUL.NE.0)             GO TO 94
#endif
      IDN = IQ(KQS+LNEXT-5) + 1
      GO TO 39

   36 IF (NQBIA.EQ.0)              GO TO 37
      LSAME = 0
      IDN   = -NQBIA
      GO TO 39

   37 IDN = IQ(KQS+LSAME-5) + 1
      GO TO 39

   38 IF (LNEXT.NE.0)  IDN=IQ(KQS+LNEXT-5)+1
   39 CONTINUE

C----              Ready I/O characteristic

      IF (ICHORG.LT.0)             GO TO 47

C--                Immediate

      IF (ICHORG.LT.8)  THEN
          NQNIO = 0
          NQIOCH(1) = ISHFTL (ICHORG, 16)
          GO TO 49
        ENDIF

C--                Copy characteristic from sister bank

      IF (ICHORG-11)         45, 43, 47
   43 IF (LSAME.EQ.0)              GO TO 45
#if !defined(CERNLIB_QDEBUG)
      NQNIO = JBYT (IQ(KQS+LSAME),19,4)
      IQLN  = LSAME - (NQNIO+IQ(KQS+LSAME-3)+1)
#endif
#if defined(CERNLIB_QDEBUG)
      NQNIO = IQNIO
#endif
      IF (NQNIO.EQ.0)  THEN
          NQIOCH(1) = LQ(KQS+IQLN)
          GO TO 49
        ELSE
          CALL UCOPY (LQ(KQS+IQLN),NQIOCH,NQNIO+1)
          NQIOSV(1) = 0
          GO TO 49
        ENDIF

C--                Find index to characteristic according to IDH

   45 LID  = LQFORM
      IF (LID.EQ.0)                GO TO 95
      LIOD = LQ(KQSP+LID-2)
      IF (NQID.LT.0)  LID=LQ(KQSP+LID)

C--                Same as last

      IF (NQID.EQ.IQ(KQSP+LID+3))  THEN
          IXIO = IQ(KQSP+LID+2)
        ELSE

C--                Search

          N = IQ(KQSP+LID+1)
          IF (N.EQ.0)              GO TO 95
          J = IUCOMP (NQID,IQ(KQSP+LID+4),N)
          IF (J.EQ.0)              GO TO 95

          LIX  = LQ(KQSP+LID-1)
          IXIO = IQ(KQSP+LIX+J)
          IQ(KQSP+LID+2) = IXIO
          IQ(KQSP+LID+3) = NQID
        ENDIF

      NQNIO = JBYT (IQ(KQSP+LIOD+IXIO+1),7,5) - 1
      GO TO 48

C--                Parameter is characteristic or index

   47 J     = JBYT (ICHORG,1,6)
      NQNIO = JBYT (ICHORG,7,5) - 1
      IOTH  = JBYT (ICHORG,12,5)
      IF (J.EQ.1)  THEN
          IF (NQNIO.NE.IOTH)       GO TO 96
          GO TO 49
        ENDIF

C--                Index

      IF (J.NE.2)                  GO TO 96
      IF (IOTH.NE.0)               GO TO 96
      IXIO = JBYT (ICHORG,17,16)
      IF (IXIO.EQ.0)               GO TO 96

      LIOD = LQ(KQSP+LQFORM-2)
      IF (IXIO.GE.IQ(KQSP+LIOD+1))    GO TO 96

C--                Same characteristic as previously lifted bank

   48 IF (IXIO.EQ.NQIOSV(1))  THEN
          NQIOCH(1) = NQIOSV(2)
          GO TO 49
        ENDIF

C--                Copy characteristic

      NQIOSV(1) = 0
      IF (NQNIO.GE.16)             GO TO 96
      CALL UCOPY (IQ(KQSP+LIOD+IXIO+1),NQIOCH,NQNIO+1)
      IOTH = JBYT (NQIOCH(1),12,5)
      IF (NQNIO.NE.IOTH)           GO TO 96
      NQIOSV(1) = IXIO
      NQIOSV(2) = NQIOCH(1)

   49 NTOT = NTOT + NQNIO

C------            Select division

      IF (JQDIVI.NE.0)             GO TO 59
      IF (LS.LT.LQSTA(KQT+1))         GO TO 58
      IF (LS.GE.LQEND(KQT+20))        GO TO 58
      IF (LS.GE.LQEND(KQT+JQDVLL))    GO TO 54
      IF (LS.LT.LQEND(KQT+2))         GO TO 57
      JQDIVI = 3
      GO TO 55

   54 JQDIVI = JQDVSY
   55 IF (LS.LT.LQEND(KQT+JQDIVI))    GO TO 61
      JQDIVI = JQDIVI + 1
      GO TO 55

   57 JQDIVI = 1
      IF (LS.LT.LQSTA(KQT+2))         GO TO 61
   58 JQDIVI = 2
      GO TO 61

   59 IF (LSAME.EQ.0)               GO TO 61
      IF (LSAME.LT.LQSTA(KQT+JQDIVI))  GO TO 97
      IF (LSAME.GE.LQEND(KQT+JQDIVI))  GO TO 97

C------            Allocate space

   61 CALL MZRESV
      NQRESV = NQRESV - NTOT
      IF (NQRESV.LT.0)             GO TO 81
      IF (JQMODE.NE.0)             GO TO 63

C--                Division is mode forward

      NQLN  = LQEND(KQT+JQDIVI)
      LE    = NQLN + NTOT
      LQEND(KQT+JQDIVI) = LE
      GO TO 65

C--                Division is mode reverse

   63 LE    = LQSTA(KQT+JQDIVI)
      NQLN  = LE - NTOT
      LQSTA(KQT+JQDIVI) = NQLN

C----              Construct bank

   65 NZ = MIN (NZERO,NQND)
      IF (NZ.EQ.0)  NZ=NQND
      IF (NZ.LT.0)  NZ=0

      NST  = NQNIO + NQNL
      NQLS = NQLN + NST + 1
      CALL VZERO (LQ(KQS+NQLN+NQNIO+1),NQNL+NZ+9)

      NQIOCH(1) = MSBYT (NST+12,NQIOCH(1),1,16)
      DO 67  J=0,NQNIO
   67 LQ(KQS+NQLN+J) = NQIOCH(J+1)

      IQ(KQS+NQLS-5) = IDN
      IQ(KQS+NQLS-4) = NQID
      IQ(KQS+NQLS-3) = NQNL
      IQ(KQS+NQLS-2) = NQNS
      IQ(KQS+NQLS-1) = NQND
      IQ(KQS+NQLS)   = ISHFTL (NQNIO,18)

C------            Set up how to link

      IF   (NQBIA-1)         72, 73, 79

C--                JBIAS -ve, insert into tree

   72 LUP   = LQSUP
      KADR  = LQSUP + NQBIA
      LNEXT = LQ(KQS+KADR)
      IF (NQBIA.NE.0)              GO TO 77
      LUP   = LQ(KQS+LUP+1)
      GO TO 77

C--                JBIAS = +1, add to linear structure

   73 LNEXT = LQSUP
      IF (LNEXT.NE.0)              GO TO 74
      LUP  = 0
      KADR = LOCF (LSUPP(1)) - LQSTOR
#if defined(CERNLIB_APOLLO)
      KADR = RSHFT (IADDR(LSUPP(1)),2) - LQSTOR
#endif
      IF (KADR.LT.LQSTA(KQT+1))       GO TO 78
      IF (KADR.LT.LQSTA(KQT+21))      GO TO 98
      GO TO 78

   74 LUP  = LQ(KQS+LNEXT+1)
      KADR = LQ(KQS+LNEXT+2)

C----              Link bank into structure

   77 IF (LNEXT.EQ.0)              GO TO 78
      LQ(KQS+NQLS)    = LNEXT
      LQ(KQS+LNEXT+2) = NQLS

   78 LQ(KQS+NQLS+1) = LUP
      LQ(KQS+NQLS+2) = KADR

      LQ(KQS+KADR)   = NQLS

   79 LP(1) = NQLS
#if defined(CERNLIB_QDEBPRI)
      IF (NQLOGL.GE.2)
     + WRITE (IQLOG,9079) JQSTOR,JQDIVI,NQLS,LQSUP,NQBIA,
     +                                  NQID,NQNL,NQNS,NQND
 9079 FORMAT (' MZLIFT-  Store/Div',2I3,' L/LSUP/JBIAS=',2I9,I6,
     F' ID,NL,NS,ND= ',A4,2I6,I9)
#endif
#include "zebra/qtrace99.inc"
      RETURN

C----              Garbage collection

   81 LQMST(KQT+1) = LQSUP
      CALL MZGAR1
      LQSUP = LQMST(KQT+1)
      IF (NQBIA.GE.1)              GO TO 61
      KADR = LOCF (LSUPP(1)) - LQSTOR
#if defined(CERNLIB_APOLLO)
      KADR = RSHFT (IADDR(LSUPP(1)),2) - LQSTOR
#endif
      IF (KADR.LT.LQSTA(KQT+1))       GO TO 83
      IF (KADR.LT.LQSTA(KQT+21))      GO TO 61
   83 LSUPP(1) = LQSUP
      GO TO 61

C----              Error conditions

   98 NQCASE = 8
      NQFATA = 1
      IQUEST(18) = KADR
      GO TO 90

   97 NQCASE = 7
      NQFATA = 1
      IQUEST(18) = LSAME
      GO TO 90

   94 NQCASE = 4
      NQFATA = 1
      IQUEST(18) = LNEXT
      GO TO 90

   96 NQCASE = 1
   95 NQCASE = NQCASE + 2
   93 NQCASE = NQCASE + 1
   92 NQCASE = NQCASE + 1
   91 NQCASE = NQCASE + 1
   90 NQFATA = NQFATA + 7
      IQUEST(11) = LQSUP
      IQUEST(12) = NQBIA
      IQUEST(13) = NQID
      IQUEST(14) = NQNL
      IQUEST(15) = NQNS
      IQUEST(16) = NQND
      IQUEST(17) = ICHORG
#include "zebra/qtofatal.inc"
      END
*      ==================================================
#include "zebra/qcardl.inc"
