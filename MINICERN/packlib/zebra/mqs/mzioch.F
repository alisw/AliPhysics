*
* $Id$
*
* $Log$
* Revision 1.2  1996/04/18 16:13:03  mclareni
* Incorporate changes from J.Zoll for version 3.77
*
* Revision 1.1.1.1  1996/03/06 10:47:21  mclareni
* Zebra
*
*
#include "zebra/pilot.h"
      SUBROUTINE MZIOCH (IODVEC,NWIOMP,CHFORM)

C-    Analyse I/O descriptor, user called, also from MZFORM / MZIOBK
C-    convert the descriptor CHFORM into the I/O characteristic
C-    stored into IODVEC of length NWIOMX

#include "zebra/zkrakc.inc"
#include "zebra/zstate.inc"
#include "zebra/zunit.inc"
#include "zebra/quest.inc"
#include "zebra/mzca.inc"
#include "zebra/zbcd.inc"
C--------------    END CDE                             --------------
      DIMENSION    IODVEC(99), NWIOMP(9)
      CHARACTER    CHFORM*(*)

      EQUIVALENCE (NGR,IQUEST(1)), (NGRU,IQUEST(2))

      DIMENSION    MU(99),             MCE(99)
      EQUIVALENCE (MU(1),IQHOLK(1)),  (MCE(1),IQCETK(1))

      DIMENSION    NBITVA(4), NBITVB(4), NBITVC(7)
      DIMENSION    MXVALA(4), MXVALB(4), MXVALC(7)

      DIMENSION    ITAB(48), INV(10)

#if defined(CERNLIB_QMVDS)
      SAVE         NBITVA, NBITVB, NBITVC
      SAVE         MXVALA, MXVALB, MXVALC
      SAVE         ITAB, INV
#endif
#if (defined(CERNLIB_QTRHOLL))&&(!defined(CERNLIB_A6M))
      DIMENSION    NAMESR(2)
      DATA  NAMESR / 4HMZIO, 4HCH   /
#endif
#if (defined(CERNLIB_QTRHOLL))&&(defined(CERNLIB_A6M))
      DATA  NAMESR / 6HMZIOCH /
#endif
#if !defined(CERNLIB_QTRHOLL)
      CHARACTER    NAMESR*8
      PARAMETER   (NAMESR = 'MZIOCH  ')
#endif

C-                 A   B   C   D   E   F   G   H   I   J   K   L   M
C-                 N   O   P   Q   R   S   T   U   V   W   X   Y   Z
C-                 0   1   2   3   4   5   6   7   8   9   +   -   *
C-                 /   (   )   $   =   b   ,   .
      DATA ITAB /       47
     +,           -1, 12, -1, 15, -1, 14, -1, 16, 13, -1, -1, -1, -1
     +,           -1, -1, -1, -1, -1, 18, -1, -1, -1, -1, -1, -1, -1
     +,            0,  1,  2,  3,  4,  5,  6,  7,  8,  9, -1, 11, 10
     +,           19, -1, -1, -1, -1, -2, -2, -2 /

      DATA INV  / 39, 38,  2,  9,  6,  4,  8, 24, 19, 40 /
C-                 *   -   B   I   F   D   H   X   S   /

      DATA  NBITVA / 32,    16,   10,   8 /
      DATA  NBITVB / 29,    14,    9,   7 /
      DATA  NBITVC / 26,    11,    6,   4, 2, 1, 1 /
      DATA  MXVALA /  0, 65536, 1024, 256 /
      DATA  MXVALB /  0, 16384,  512, 128 /
      DATA  MXVALC /  0,  2048,   64,  16, 4, 2, 2 /

#include "zebra/q_sbit1.inc"
#include "zebra/q_sbyt.inc"

#include "zebra/qtrace.inc"

      NWIOMX = NWIOMP(1)
      NCH    = LEN (CHFORM)
      IF (NCH.GE.121)              GO TO 90
      CALL UCTOH1 (CHFORM,IQHOLK,NCH)

      CALL IZBCDT (NCH,ITAB)
      NCH = IQUEST(1)
      IF (IQUEST(2).NE.0)          GO TO 91
      IF (IQUEST(1).EQ.0)          GO TO 91

C--------          SCAN CHARACTER BY CHARACTER

C-       NUMERIC   0 ... 9  10 11  12 13 14 15 16 17 18 19
C-           FOR   0 ... 9   *  -   B  I  F  D  H     S  /

      JPOSR  = -1
      JPOSIN = -1
      IVAL = 0
      JCH  = 0
      JU   = 0
   21 NVAL = 0
   22 JCH  = JCH + 1
      NUM  = MCE(JCH)
      IF (NUM.GE.10)               GO TO 24

C--                NUMERIC

      NVAL = 10*NVAL + NUM
      IF (JCH.LT.NCH)              GO TO 22
      GO TO 92

C--                COUNT =  * OR -

   24 IF (NUM.GE.12)               GO TO 26
      IF (NVAL.NE.0)               GO TO 92
      IF (NUM.EQ.11)  THEN
          NVAL   = -1
          JPOSIN = JU
          IF (JPOSR.GE.0)          GO TO 93
        ENDIF
      IF (JCH.EQ.NCH)              GO TO 92
      JCH = JCH + 1
      NUM = MCE(JCH)
      IF (NUM.LT.12)               GO TO 92
      IF (NUM.GE.19)               GO TO 92
      GO TO 27

C--                TYPE LETTER

   26 IF (NUM.EQ.19)               GO TO 29
      IF (NUM.EQ.18)               GO TO 92
      IF (NVAL.EQ.0)               GO TO 92
      IF (NUM.EQ.15)  THEN
          IF (NVAL.NE.2*(NVAL/2))    GO TO 92
        ENDIF
      IVAL = 7
   27 MU(JU+1) = NUM - 11
      MU(JU+2) = NVAL
      JU = JU + 2
      IF (JCH.EQ.NCH)              GO TO 31
      IF (JPOSIN.LT.0)             GO TO 21
      GO TO 94

C--                REPEAT MARK

   29 IF (NVAL.NE.0)               GO TO 92
      IF (JPOSR.GE.0)              GO TO 95
      IF (JCH.EQ.NCH)              GO TO 92
      JPOSR = JU
      GO TO 21

C--------          NO TRAILING REGIONS

   31 NU    = JU
      NSECA = NU/2
      JU    = 2
      IOWD  = 65
      NWIO  = 0
      JFL12 = 1
      IF (JPOSR.GE.0)  THEN
          IF (JPOSR+2.NE.NU)       GO TO 41
          IF (MU(NU-1).NE.7)       GO TO 41
C-        continue if  'CT ... CT / *S'
        ELSE
          IF (MU(NU).EQ.0)  JFL12=2
        ENDIF

C--        FORMATS TO BE HANDLED IN CLASSES 1 AND 2 (OR 0)
C-    JFL12 = 2 :  '... *T'
C-            1 :  '... -T'
C-                 '... NT'    AS  '... -T'
C-                 '... / *S'  AS  '... -S'
C-                 '/ NT'      AS  '-T'

   32 NSECA = NSECA - 1

C----              CLASS 0 : '-T'   OR   '*T'

      IQUEST(12) = MU(NU-1)
      IF   (NSECA.EQ.0)   THEN
          IF (JFL12.EQ.1)          GO TO 82
          IQUEST(12) = MSBIT1 (IQUEST(12),4)
          GO TO 82
        ENDIF

C----              CLASS 0 : 'CT -T'   OR   'CT *T'

      IF (NSECA.GE.2)              GO TO 33
      IF (MU(2).GE.64)             GO TO 34
      IF (JFL12.EQ.2)  IQUEST(12)= MSBIT1(IQUEST(12),4)
      IQUEST(12) = MSBYT (MU(1),IQUEST(12),5,3)
      IQUEST(12) = MSBYT (MU(2),IQUEST(12),8,6)
      GO TO 82

C----              CLASS 1 : 'CT CT CT -T'
C-                 CLASS 2 : 'CT CT CT *T'

   33 IF (IVAL+NSECA.EQ.2)         GO TO 38
   34 IQUEST(12) = MSBYT (JFL12,IQUEST(12),14,2)
      IQUEST(12) = MSBYT (MU(1),IQUEST(12),5,3)
      IQUEST(13) = MU(2)
      JBTF  = 8
      IF (NSECA.GE.4)              GO TO 36
      IOWD  = 2177
      NWIO  = 1
      IF (NSECA.EQ.1)              GO TO 82

      NGR = NSECA
      CALL MZIOCF (0,MXVALA)
      IF (NGR.NE.NGRU)             GO TO 36
      NBT   = NBITVA(NGRU)
      GO TO 71

C----              CLASS 1 : 'CT ... CT -T'
C-                 CLASS 2 : 'CT ... CT *T'

   36 IQUEST(12) = MSBIT1 (IQUEST(12),4)
      NGR = MIN (NSECA,3)
      CALL MZIOCF (0,MXVALB)
      NBT   = NBITVB(NGRU)
      GO TO 70

C----              CLASS 1 : '*T *T -T'
C-                 CLASS 2 : '*T *T *T'

   38 IQUEST(12) = 16*IQUEST(12)
      IQUEST(12) = MSBYT (MU(1),IQUEST(12), 8,3)
      IQUEST(12) = MSBYT (MU(3),IQUEST(12),11,3)
      IQUEST(12) = MSBYT (JFL12,IQUEST(12),14,2)
      GO TO 82

C--------          WITH TRAILING REGIONS

   41 NSECL = JPOSR/2
      IF (NSECA.GE.3)              GO TO 44

C----              CLASS 3 : '/ *T'

      IQUEST(12) = MU(NU-1)
      IQUEST(12) = MSBYT (3,IQUEST(12),14,2)
      IF (NSECA.EQ.2)              GO TO 42
      IF (MU(2).EQ.0)              GO TO 82
      GO TO 32

C----              CLASS 3 : 'CT / *T'   OR   '/ CT *T'

   42 IF (MU(4).NE.0)              GO TO 44
      IF (MU(2).GE.64)             GO TO 44
      IQUEST(12) = MSBYT (MU(1),IQUEST(12),5,3)
      IQUEST(12) = MSBYT (MU(2),IQUEST(12),8,6)
      IF (NSECL.EQ.1)              GO TO 82
      IQUEST(12) = MSBIT1 (IQUEST(12),4)
      GO TO 82

C----              CLASS 4 :  'CT / CT CT CT'  OR  'CT CT / CT CT'

   44 IF (NSECL.EQ.0)              GO TO 51
      IF (NSECL.GE.3)              GO TO 61
      IF (NSECA.GE.5)              GO TO 61
      IF (IVAL+NSECA.EQ.3)         GO TO 48

      NGR = NSECA
      CALL MZIOCF (0,MXVALA)
      IF (NGR.NE.NGRU)             GO TO 61

      IQUEST(12) = MU(1)
      IQUEST(13) = MU(2)
      IF (NSECL.EQ.2)  IQUEST(12)=IQUEST(12)+8
      IQUEST(12) = MSBIT1 (IQUEST(12),16)

      JBTF  = 5
      NBT   = NBITVA(NGRU)
      IOWD  = 2177
      NWIO  = 1
      GO TO 71

C----              CLASS 4 :  '*T / *T *T'   OR   '*T *T / *T'

   48 IQUEST(12) = 8*(2*MU(1)+NSECL-1)
      IQUEST(12) = MSBYT (MU(3),IQUEST(12), 8,3)
      IQUEST(12) = MSBYT (MU(5),IQUEST(12),11,3)
      IQUEST(12) = MSBIT1 (IQUEST(12),16)
      GO TO 82

C--------          NO LEADING REGIONS

C----              CLASS 5 : '/ CT CT CT CT'

   51 IF (IVAL+NSECA.EQ.3)         GO TO 58
      IQUEST(12) = MU(1)
      IQUEST(13) = MU(2)
      IQUEST(12) = MSBYT (5,IQUEST(12),14,3)
      JBTF  = 5
      IF (NSECA.GE.5)              GO TO 55

      NGR = NSECA
      CALL MZIOCF (0,MXVALA)
      IF (NGR.NE.NGRU)             GO TO 55
      NBT   = NBITVA(NGRU)
      IOWD  = 2177
      NWIO  = 1
      GO TO 71

C----              CLASS 5 : '/ CT ... CT'

   55 IQUEST(12) = MSBIT1 (IQUEST(12),4)
      NGR = MIN (NSECA,4)
      CALL MZIOCF (0,MXVALB)
      NBT  = NBITVB(NGRU)
      GO TO 70

C----              CLASS 5 : '/ *T *T *T'

   58 IQUEST(12) = 16*MU(1)
      IQUEST(12) = MSBYT (MU(3),IQUEST(12), 8,3)
      IQUEST(12) = MSBYT (MU(5),IQUEST(12),11,3)
      IQUEST(12) = MSBYT     (5,IQUEST(12),14,3)
      GO TO 82

C--------          CLASS 6 : 'CT ... CT / CT ... CT'

   61 IQUEST(12) = NSECL
      IQUEST(13) = MU(2)
      IF (NSECL.GE.16)             GO TO 96
      IQUEST(12) = MSBYT (MU(1),IQUEST(12),5,3)
      IQUEST(12) = MSBYT (6,IQUEST(12),14,3)

      JBTF  = 8
      NGR   = 3
      CALL MZIOCF (0,MXVALB)
      NBT   = NBITVB(NGRU)

C------            LONG PACKING

C--                PACK FIRST I/O WORD

   70 IF (NGRU.EQ.1)               GO TO 73

   71 JBTC = 1
      DO  72  JL=2,NGRU
      IQUEST(12) = MSBYT (MU(JU+1),IQUEST(12),JBTF,3)
      JBTF  = JBTF + 3
      JBTC  = JBTC + NBT
      IQUEST(13) = MSBYT (MU(JU+2),IQUEST(13),JBTC,NBT)
   72 JU = JU + 2

      IF (NGRU.EQ.NSECA)           GO TO 82
   73 NSECD = NGRU
      JWIO  = 13

C----              PACK NEXT I/O WORD

   74 JWIO = JWIO + 1
      IQUEST(JWIO) = MU(JU+1)
      JBT  = 4
      NGRU = 1
      NGR  = MIN (7,NSECA-NSECD)
      IF (NGR.EQ.1)                GO TO 77

      CALL MZIOCF (JU,MXVALC)
      IF (NGRU.EQ.1)               GO TO 77
      JUST = JU
      DO  76  JL=2,NGRU
      JU  = JU + 2
      IQUEST(JWIO) = MSBYT (MU(JU+1),IQUEST(JWIO),JBT,3)
   76 JBT = JBT + 3
      JU  = JUST

   77 IQUEST(JWIO-1) = MSBYT (NGRU,IQUEST(JWIO-1),30,3)
      NBT  = NBITVC(NGRU)

      DO  79  JL=1,NGRU
      IQUEST(JWIO) = MSBYT (MU(JU+2),IQUEST(JWIO),JBT,NBT)
      JBT  = JBT + NBT
   79 JU  = JU  + 2

      NSECD = NSECD + NGRU
      IF (NSECD.LT.NSECA)          GO TO 74
      NWIO = JWIO - 12

      IF (NWIO.GE.NWIOMX)          GO TO 97
      IF (NWIO.GE.16)              GO TO 97
      IOWD  = 64*(32*NWIO+NWIO+1) + 1

   82 IOWD = MSBYT (IQUEST(12),IOWD,17,16)
      IQUEST(12) = IOWD
      IQUEST(1)  = NWIO

      CALL UCOPY (IQUEST(12),IODVEC,NWIO+1)
      IQCETK(121) = IQBLAN
#if defined(CERNLIB_QDEBPRI)
      IF (NQLOGM.GE.1)  WRITE (IQLOG,9088) NWIO,CHFORM
 9088 FORMAT (' MZIOCH-',I5,' extra I/O words for Format ',A)
#endif
#include "zebra/qtrace99.inc"
      RETURN

C------            ERROR CONDITIONS

C--                CHARACTER STRING TOO LONG
   90 NQFATA = 2
      IQUEST(12) = NCH
      GO TO 99

C--                STRING INVALID
   91 NQCASE = 1
      NQFATA = 3
      IQUEST(12) = IQUEST(1)
      IQUEST(13) = IQUEST(2)
      IF (IQUEST(1).EQ.0)          GO TO 99
      GO TO 98

C--                I/O DESCRIPTOR TOO LONG
   97 NQCASE = 7
      IQUEST(12) = NWIOMX
      IQUEST(13) = NWIO + 1
      GO TO 98

C--                TOO MANY SECTORS
   96 NQCASE = 6
      IQUEST(12) = NSECA
      IQUEST(13) = NSECL
      GO TO 98

C--                / OCCURS TWICE
   95 NQCASE = 1

C--                -T NOT LAST SECTOR IN THE STRING
   94 NQCASE = NQCASE + 1

C--                -T NOT ALLOWED IN REPEAT
   93 NQCASE = NQCASE + 1

C--                BAD SYNTAX
   92 NQCASE = NQCASE + 2
      IQUEST(12) = JCH
      IQUEST(13) = 0

   98 DO  88  JCH=1,NCH
      JCET = MCE(JCH)
      IF (JCET.LT.10) THEN
          MCE(JCH)=IQNUM(JCET+1)
        ELSE
          JCET = INV(JCET-9)
          MCE(JCH) = IQLETT(JCET)
        ENDIF
   88 CONTINUE
      CALL UTRANS (MCE,IQUEST(14),NCH,1,4)
      NQFATA = (NCH-1)/4 + 4

   99 IQUEST(11) = IQCETK(121)
#include "zebra/qtofatal.inc"
      END
*      ==================================================
#include "zebra/qcardl.inc"
