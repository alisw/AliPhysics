*
* $Id$
*
* $Log$
* Revision 1.2  1996/04/18 16:10:24  mclareni
* Incorporate changes from J.Zoll for version 3.77
*
* Revision 1.1.1.1  1996/03/06 10:47:12  mclareni
* Zebra
*
*
#include "zebra/pilot.h"
#if defined(CERNLIB_FZALFA)
      SUBROUTINE FZIASC (NRSKPP)

C-    Read one physical record from cards in ALFA mode,
C-    service routine to FZIN, called via FZIPHA

C-    Do not expand fast blocks if NRSKIP not zero

C-    The record is stored into words LQ(LBUFA) to LQ(LBUFE-1)
C-    It is a 'fast' record if JFAST .NE. 0

#include "zebra/zbcd.inc"
#include "zebra/zkrakc.inc"
#include "zebra/zstate.inc"
#include "zebra/zunit.inc"
#include "zebra/mqsysh.inc"
#include "zebra/fzci.inc"
C--------------    End CDE                             --------------
      DIMENSION    NRSKPP(9)
      DIMENSION    INITV(6), NBV(6), ICHSUM(2,2)

      EQUIVALENCE   (LUN,IQUEST(90)), (JSKIP,IQUEST(91))
      EQUIVALENCE (MRSTA,IQUEST(92)), (MREND,IQUEST(93))
      EQUIVALENCE  (JTKC,IQUEST(96)),  (JTKL,IQUEST(97))
     +,            (JTKE,IQUEST(98))

#if defined(CERNLIB_QMVDS)
      SAVE         INITV
#endif
#if (defined(CERNLIB_QTRHOLL))&&(!defined(CERNLIB_A6M))
      DIMENSION    NAMESR(2)
      DATA  NAMESR / 4HFZIA, 4HSC   /
#endif
#if (defined(CERNLIB_QTRHOLL))&&(defined(CERNLIB_A6M))
      DATA  NAMESR / 6HFZIASC  /
#endif
#if !defined(CERNLIB_QTRHOLL)
      CHARACTER    NAMESR*8
      PARAMETER   (NAMESR = 'FZIASC  ')
#endif
      DATA  INITV / 0, 1, 2, 3, 0, 134217727 /
C-                                 = 7FFFFFF Hex

#include "zebra/q_or.inc"
#include "zebra/q_shiftl.inc"
#include "zebra/q_jbyt.inc"

#include "zebra/qtraceq.inc"

      LUN    = LUNI
      LBUFA  = L4STAI
      LBUFE  = L4STAI + MAXREI
      NRSKIP = NRSKPP(1)

      LBUFC  = LBUFA
      IFLREP = 0
      JCHSUM = 0

      ICHSUM(1,1) = 0
      ICHSUM(2,1) = 0

C----              Read first line of next physical record

      IQUEST(1) = 0
      JTKC  = 81
      JTKE  = JTKC
      MRSTA = 0
      MREND = 0
      JSKIP = 7

      CALL FZIALN
      IF (IQUEST(1).NE.0)          GO TO 101
      IF (MRSTA.GT.1)  THEN
          JFAST = 0
          LBUFC = LBUFC+4
        ELSE
          JFAST = 7
          IF (NRSKIP.NE.0)         GO TO 999
          IF (NFASTI.EQ.0)         GO TO 999
        ENDIF

      MRSTA = 0
      JSKIP = 0
      ITYPE = 0

C-------           Do next word

   31 IF (JTKC.GE.JTKL)  THEN
          IF (MREND.NE.0)          GO TO 911
          CALL FZIALN
          IF (IQUEST(1).NE.0)      GO TO 101
        ENDIF

C--                Handle running type, check termination

      IF (ITYPE.LE.0)              GO TO 34
      IF (IQKRAK(JTKC).LT.33)      GO TO 41
      GO TO 70

C--                Find type of next word

   34 JTYPC = IQKRAK(JTKC)
      JTKC  = JTKC + 1
#if defined(CERNLIB_QDEVZE)
      IF (LOGLVI.GE.6)  WRITE (IQLOG,9834) JTYPC,IQLETT(JTYPC)
 9834 FORMAT (' FZIASC-  Next control code/char.=',I4,1X,A1)
#endif
      IF (JTYPC.GE.42)             GO TO 71
      IF   (JTYPC-26)        36, 912, 35
   35 IWORD = JTYPC - 27
      IF (IWORD.LT.10)             GO TO 61
      JTYPC = JTYPC - 11
   36 JTYPC = JTYPC - 1
      JTYPE = JTYPC / 5
      JTYPS = JTYPC - 5*JTYPE
      JTYPE = JTYPE + 1

      NBUSE = 6 - JTYPS
      IF (JTYPE.GE.5)  NBUSE=NBUSE-1

      INIT = ISHFTL (INITV(JTYPE), 5)

C----              Compose next word

C--                copy bytes and check validity

   41 DO 42  J=1,NBUSE
      NBV(J) = IQKRAK(JTKC) - 1
      IF (NBV(J).GE.32)            GO TO 913
   42 JTKC = JTKC + 1

      IWORD = IOR (NBV(1),INIT)

      IF (JTYPE.LT.5)              GO TO 51

C--                Compose NBUSE 5-bit bytes for JTYPE=5,6

      GO TO ( 60, 44, 45, 46, 47, 48), NBUSE

   44 IWORD = IOR (NBV(2), ISHFTL(IWORD,5))
      GO TO 60

   45 IWORD = IOR (NBV(2), ISHFTL(IWORD,5))
      IWORD = IOR (NBV(3), ISHFTL(IWORD,5))
      GO TO 60

   46 IWORD = IOR (NBV(2), ISHFTL(IWORD,5))
      IWORD = IOR (NBV(3), ISHFTL(IWORD,5))
      IWORD = IOR (NBV(4), ISHFTL(IWORD,5))
      GO TO 60

   47 IWORD = IOR (NBV(2), ISHFTL(IWORD,5))
      IWORD = IOR (NBV(3), ISHFTL(IWORD,5))
      IWORD = IOR (NBV(4), ISHFTL(IWORD,5))
      IWORD = IOR (NBV(5), ISHFTL(IWORD,5))
      GO TO 60

   48 IWORD = IOR (NBV(2), ISHFTL(IWORD,5))
      IWORD = IOR (NBV(3), ISHFTL(IWORD,5))
      IWORD = IOR (NBV(4), ISHFTL(IWORD,5))
      IWORD = IOR (NBV(5), ISHFTL(IWORD,5))
      IWORD = IOR (NBV(6), ISHFTL(IWORD,5))
      GO TO 60

C--                Compose NBUSE 5-bit bytes for JTYPE=1,2,3,4

   51 GO TO ( 60, 54, 55, 56, 57, 58), NBUSE

   54 IWORD = IOR (NBV(2), ISHFTL(IWORD,5))
      IWORD = ISHFTL (IWORD,20)
      GO TO 60

   55 IWORD = IOR (NBV(2), ISHFTL(IWORD,5))
      IWORD = IOR (NBV(3), ISHFTL(IWORD,5))
      IWORD = ISHFTL (IWORD,15)
      GO TO 60

   56 IWORD = IOR (NBV(2), ISHFTL(IWORD,5))
      IWORD = IOR (NBV(3), ISHFTL(IWORD,5))
      IWORD = IOR (NBV(4), ISHFTL(IWORD,5))
      IWORD = ISHFTL (IWORD,10)
      GO TO 60

   57 IWORD = IOR (NBV(2), ISHFTL(IWORD,5))
      IWORD = IOR (NBV(3), ISHFTL(IWORD,5))
      IWORD = IOR (NBV(4), ISHFTL(IWORD,5))
      IWORD = IOR (NBV(5), ISHFTL(IWORD,5))
      IWORD = ISHFTL (IWORD,5)
      GO TO 60

   58 IWORD = IOR (NBV(2), ISHFTL(IWORD,5))
      IWORD = IOR (NBV(3), ISHFTL(IWORD,5))
      IWORD = IOR (NBV(4), ISHFTL(IWORD,5))
      IWORD = IOR (NBV(5), ISHFTL(IWORD,5))
      IWORD = IOR (NBV(6), ISHFTL(IWORD,5))

C--                Store composed word

#if defined(CERNLIB_QDEVZE)
   60 IF (LOGLVI.GE.6)  WRITE (IQLOG,9860) IWORD,(NBV(J),J=1,NBUSE)
 9860 FORMAT (' FZIASC-  Composed word= ',Z8,6I3)
      GO TO 62
   61 IF (LOGLVI.GE.6)  WRITE (IQLOG,9861) IWORD
 9861 FORMAT (' FZIASC-   Integer word= ',I8)
   62 CONTINUE
#endif
#if !defined(CERNLIB_QDEVZE)
   60 CONTINUE
   61 CONTINUE
#endif
      IF (JCHSUM.NE.0)             GO TO 82

      ICHSUM(1,1) = ICHSUM(1,1) + JBYT (IWORD,17,16)
      ICHSUM(2,1) = ICHSUM(2,1) + JBYT (IWORD, 1,16)

      IF (IFLREP.NE.0)             GO TO 63
      IF (LBUFC.GE.LBUFE)          GO TO 914
      LQ(LBUFC) = IWORD
      LBUFC = LBUFC + 1
      IF (ITYPE.GE.0)              GO TO 31
      ITYPE = JTYPE
      GO TO 31

C----              Repetition executed

   63 N     = IWORD + 1
      LBUFN = LBUFC + N
      IF (LBUFN.GT.LBUFE)          GO TO 915
      IF (N.LE.0)                  GO TO 916

      DO 64  L=LBUFC,LBUFN-1
   64 LQ(L) = IWDREP
      LBUFC = LBUFN

      IFLREP = 0
      GO TO 31

C----              Control symbols

C--                = : repetition signalled

   70 JTYPC = IQKRAK(JTKC)
      JTKC  = JTKC + 1
   71 IF   (JTYPC-44)        912, 72, 74
   72 IF (JCHSUM.EQ.1)             GO TO 920
      IWDREP = IWORD
      IFLREP = 7
      ITYPE  = 0
      GO TO 31

C--                close sq bracket : stop running type

   74 IF   (JTYPC-54)        912, 75, 77
   75 ITYPE = 0
      GO TO 31

C--                open sq bracket : start running type

   77 IF   (JTYPC-58)        912, 78, 81
   78 IF (ITYPE.LT.0)              GO TO 917
      ITYPE = -7
      GO TO 31

C--                < : end of physical record

   81 IF (JTYPC.NE.60)             GO TO 912
      IF (LBUFC.NE.LBUFE)          GO TO 918
      IF (JCHSUM.NE.0)             GO TO 999
      IF (IFLREP.NE.0)             GO TO 920
      ITYPE  = 0
      JCHSUM = 1
      GO TO 31

C--                Check-sum reading

   82 IF (IFLREP.NE.0)  IWORD = IWDREP
      ICHSUM(JCHSUM,2) = IWORD

      JCHSUM = JCHSUM + 1
      IF (JCHSUM.EQ.2)             GO TO 31

      IF (ICHSUM(1,1).NE.ICHSUM(1,2))    GO TO 919
      IF (ICHSUM(2,1).NE.ICHSUM(2,2))    GO TO 919

#include "zebra/qtrace99.inc"
      IQUEST(92) = JFAST
      RETURN

C------            Errors

  101 IF (IQUEST(1).LT.0)          GO TO 999
      IF (IQUEST(1).EQ.7799)       GO TO 910

C--                Read Error

      JRETCD = 7
      JERROR = 309
      NWERR  = 1
      IQUEST(14) = IQUEST(1)
      GO TO 999

C----              Record context errors

C-    JERROR = 310  Invalid character in column 1
  910 JERROR = -1

C-    JERROR = 311  Record shorter than expected
  911 JERROR = JERROR - 1

C-    JERROR = 312  Faulty type code
  912 JERROR = JERROR - 1

C-    JERROR = 313  Faulty numeric value, > 31
  913 JERROR = JERROR - 1

C-    JERROR = 314  Record longer than expected
  914 JERROR = JERROR - 1

C-    JERROR = 315  Repetition count overshoots record end
  915 JERROR = JERROR - 1

C-    JERROR = 316  Repetition count negative
  916 JERROR = JERROR - 1

C-    JERROR = 317  Double open square bracket
  917 JERROR = JERROR - 1

C-    JERROR = 318  Record shorter than expected
  918 JERROR = JERROR - 1

C-    JERROR = 319  Check-sum error
  919 JERROR = JERROR - 1

C-    JERROR = 320  Illegal combination  =<  or  <=
  920 JERROR = JERROR + 320
      JRETCD = 5
      NWERR  = 0
      IQUEST(1) = 7
      GO TO 999
      END
*      ==================================================
#include "zebra/qcardl.inc"
#endif
