*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/02/15 17:49:42  mclareni
* Kernlib
*
*
#include "kerngen/pilot.h"
      SUBROUTINE CKRACK (CHV,JLP,JRP,IFLAGD)
C
C CERN PROGLIB# M432    CKRACK          .VERSION KERNFOR  4.29  910718
C ORIG. 12/06/91, JZ
C
C-    Read the next number from CHV(JL:JR)
C-    formats: (1) bits     - #On octal, or #Bn binary, or #Xn hex
C-             (2) integer  - [+|-]n
C-             (3) floating - [+|-][n][.][f][E][+|-][n]
C-             (4) double   - [+|-][n][.][f]D[+|-][n]
C-
C-    Returns: NDSLAT   number of numeric digits seen
C-             NESLAT   COL(NESLAT) is the terminating character
C-             NFSLAT   mode: -ve bad, 0 blank, 1 B, 2 I, 3 F, 4 D
C-             NGSLAT   = zero if correct termination
C-                      = NESLAT otherwise
C-             NUM(1) or ANUM(1) or DNUM  returns the value

      DIMENSION    JLP(9), JRP(9)
      CHARACTER    CHV(512)*1

      COMMON /SLATE/ NDIGT,NESLAT,MODE,NGSLAT,NUM(2)
     +,              IVALV(6),NEXPV(6), JXOP,JXME,JXFA,JXFE
     +,              JTERM,IPHASE, NEGM,NEGE,NEXPM, NERR, DUMMY(12)

      REAL             ANUM(2)
      DOUBLE PRECISION DNUM, DFRACT
      EQUIVALENCE (ANUM(1),NUM(1))
      EQUIVALENCE (DNUM,NUM(1))

      DIMENSION    ISLATE(40)
      EQUIVALENCE (ISLATE(1),NDIGT)

#include "kerngen/wordsize.inc"
#if !defined(CERNLIB_QISASTD)
#include "kerngen/q_shift.inc"
#endif

      JJ = JLP(1)
      JR = JRP(1)

      DO 12  J=1,28
   12 ISLATE(J) = 0

C----              Look at the first character of the number

   17 IF (JJ.GT.JR)                GO TO 90
      IF (CHV(JJ).EQ.' ')  THEN
          JJ = JJ + 1
          GO TO 17
        ELSEIF (CHV(JJ).EQ.'#')  THEN
          GO TO 71
        ELSEIF (CHV(JJ).EQ.'+')  THEN
          JJ = JJ + 1
        ELSEIF (CHV(JJ).EQ.'-')  THEN
          NEGM = 7
          JJ  = JJ + 1
        ENDIF

C----              Read an integer

   21 JTERM  = 0
      NDIG   = 0
      IVAL   = 0
   22 IF (JJ.GT.JR)          GO TO 27
#if defined(CERNLIB_QASCII)
      K = ICHAR (CHV(JJ))
      K = K - 48
      IF (K.LT.0)            GO TO 26
      IF (K.GE.10)           GO TO 26
#endif
#if defined(CERNLIB_QEBCDIC)
      K = ICHAR (CHV(JJ))
      K = K - 240
      IF (K.LT.0)            GO TO 26
      IF (K.GE.10)           GO TO 26
#endif
#if !defined(CERNLIB_QISASTD)
      IF (ISHFTR(IVAL,NBITPW-5).NE.0)   GO TO 24
#endif
#if defined(CERNLIB_QISASTD)
      IF (ISHFT(IVAL,5-NBITPW).NE.0)    GO TO 24
#endif
   23 JJ   = JJ + 1
      IVAL = 10*IVAL + K
      NDIG = NDIG + 1
      GO TO 22

C--       getting near the integer capacity

   24 IF (IPHASE.NE.0)             GO TO 93
      IF (JXOP.LT.6)  THEN
          JXOP = JXOP + 1
          IVALV(JXOP) = IVAL
          NEXPV(JXOP) = NDIG
        ENDIF
      IVAL = 0
      GO TO 23

   26 JTERM  = INDEX (' .+-EDed', CHV(JJ)) - 1
C-                     01234567

   27 IF (NDIG.NE.0)  THEN
          NDIGT = NDIGT + NDIG
          IF (JXOP.LT.6)  THEN
              JXOP = JXOP + 1
              IVALV(JXOP) = IVAL
              NEXPV(JXOP) = NDIG
            ENDIF
        ENDIF

C--       IPHASE = 0 :  IVAL is the leading integer
C-                 2 :       is the exponent integer

      IF (IPHASE.NE.0)             GO TO 51
      JXME = JXOP
      IF (JTERM.GE.1)              GO TO 31

C--       pure integer

   28 IF (NDIGT.EQ.0)              GO TO 91
      MODE = 2
      NUM(1) = IVALV(1)
      IF (JXME.LT.2)               GO TO 29
      N = NEXPV(2) - NEXPV(1)
      IF (N.GE.2)                  GO TO 92
#if !defined(CERNLIB_QISASTD)
      IVALV(1) = ISHFTL(NUM(1),2) + IVALV(1)
      NUM(1)   = ISHFTL(IVALV(1),1)
#endif
#if defined(CERNLIB_QISASTD)
      IVALV(1) = ISHFT(NUM(1),2) + IVALV(1)
      NUM(1)   = ISHFT(IVALV(1),1)
#endif
      IF (NUM(1).LT.0)             GO TO 92
      NUM(1) = NUM(1) + IVALV(2)
      IF (NUM(1).LT.0)             GO TO 92
   29 IF (NEGM.NE.0) NUM(1) = -NUM(1)
   30 NESLAT = JJ
      IF (JTERM.EQ.0)              RETURN
      NGSLAT = JJ
      RETURN

C----              Read the fraction

   31 IF (JTERM.NE.1)              GO TO 41
      JXFA = JXOP
      MODE = 3
      JJ   = JJ + 1

      JTERM  = 0
      NDIG   = 0
   32 IVAL   = 0
      IVALNT = 0
      NTRAIL = 0
   33 IF (JJ.GT.JR)          GO TO 37
#if defined(CERNLIB_QASCII)
      K = ICHAR (CHV(JJ))
      K = K - 48
      IF (K.LT.0)            GO TO 36
      IF (K.GE.10)           GO TO 36
#endif
#if defined(CERNLIB_QEBCDIC)
      K = ICHAR (CHV(JJ))
      K = K - 240
      IF (K.LT.0)            GO TO 36
      IF (K.GE.10)           GO TO 36
#endif
#if !defined(CERNLIB_QISASTD)
      IF (ISHFTR(IVAL,NBITPW-5).NE.0)   GO TO 34
#endif
#if defined(CERNLIB_QISASTD)
      IF (ISHFT(IVAL,5-NBITPW).NE.0)    GO TO 34
#endif
      JJ   = JJ + 1
      IVAL = 10*IVAL + K
      NDIG = NDIG + 1
      IF (K.EQ.0)  THEN
          NTRAIL = NTRAIL + 1
        ELSE
          NTRAIL = 0
          IVALNT = IVAL
        ENDIF
      GO TO 33

C--       getting near the integer capacity

   34 IF (JXOP.LT.6)  THEN
          JXOP = JXOP + 1
          IVALV(JXOP) = IVALNT
          NEXPV(JXOP) = NTRAIL - NDIG
        ENDIF
      GO TO 32

   36 JTERM  = INDEX (' .+-EDed', CHV(JJ)) - 1
C-                     01234567

   37 NDIGT = NDIGT + NDIG
      IF (IVAL.NE.0)  THEN
          IF (JXOP.LT.6)  THEN
              JXOP = JXOP + 1
              IVALV(JXOP) = IVALNT
              NEXPV(JXOP) = NTRAIL - NDIG
            ENDIF
        ENDIF
      JXFE = JXOP
      IF (JTERM.LT.2)              GO TO 52

C----              Read the exponent

   41 IPHASE = 2
      IF (NDIGT.EQ.0)              GO TO 91
      IF (JTERM.GE.4)              GO TO 44
      IF (JJ.EQ.JR)  THEN
          IF (MODE.EQ.0)           GO TO 28
          GO TO 52
        ENDIF
      MODE = 3
      NEGE = JTERM - 2
      JJ   = JJ + 1
      GO TO 21

   44 IF (JTERM.GE.6)  JTERM = JTERM - 2
      MODE = JTERM - 1
      JJ   = JJ + 1
      IF (JJ.GT.JR)  THEN
          JTERM = 0
          GO TO 52
        ENDIF

      J = INDEX ('+-', CHV(JJ))
      IF (J.EQ.0)                  GO TO 21
      IF (JJ.EQ.JR)                GO TO 52
      NEGE = J - 1
      JJ   = JJ + 1
      GO TO 21

C--                Exponent complete, construct the number

   51 NEXPM = IVAL
      IF (NEGE.NE.0)  NEXPM = -NEXPM
   52 IF (IFLAGD.GE.0)  THEN
          IF (MODE.EQ.4)           GO TO 61
          IF (IFLAGD.NE.0)         GO TO 61
        ENDIF
      ANUM(1) = 0.
      ANUM(2) = 0.
      IF (JXME.EQ.0)               GO TO 56

C--       single precision, integer part

      ANUM(1) = REAL(IVALV(1))
      IF (JXME.GE.2)  THEN
          N = NEXPV(2) - NEXPV(1)
          DO 53  J=1,N
   53     ANUM(1) = 10. * ANUM(1)
          ANUM(1) = ANUM(1) + REAL(IVALV(2))
        ENDIF

      IF (NEXPM.EQ.0)              GO TO 56
      IF (NEXPM.LT.0)              GO TO 55
      IF (NEXPM.GT.9)              GO TO 55
      DO 54  J=1,NEXPM
   54 ANUM(1) = ANUM(1) * 10.
      GO TO 56

   55 ANUM(1) = ANUM(1) * 10.**NEXPM

C--       single precision, fraction

   56 IF (JXFE.LE.JXFA)            GO TO 60
      DO 59  JXOP=JXFA+1,JXFE
      ANUM(2) = REAL(IVALV(JXOP))
      NEXPU   = NEXPV(JXOP) + NEXPM
      IF (NEXPU.EQ.0)              GO TO 59
      IF (NEXPU.LT.0)              GO TO 58
      IF (NEXPU.GT.9)              GO TO 58
      DO 57  J=1,NEXPU
   57 ANUM(2) = ANUM(2) * 10.
      GO TO 59

   58 ANUM(2) = ANUM(2) * 10.**NEXPU
   59 ANUM(1) = ANUM(1) + ANUM(2)
   60 IF (NEGM.NE.0)  ANUM(1) = -ANUM(1)
      GO TO 30

C--       double precision, integer part

   61 MODE = 4
      DNUM = 0.
      IF (JXME.EQ.0)               GO TO 66

      DNUM = DBLE(IVALV(1))
      IF (JXME.GE.2)  THEN
          N = NEXPV(2) - NEXPV(1)
          DO 63  J=1,N
   63     DNUM = DNUM * 10.
          DNUM = DNUM + DBLE(IVALV(2))
        ENDIF

      IF (NEXPM.EQ.0)              GO TO 66
      IF (NEXPM.LT.0)              GO TO 65
      IF (NEXPM.GT.9)              GO TO 65
      DO 64  J=1,NEXPM
   64 DNUM = DNUM * 10.
      GO TO 66

   65 DNUM = DNUM * 10.D0**NEXPM

C--       double precision, fraction

   66 IF (JXFE.LE.JXFA)            GO TO 70
      DO 69  JXOP=JXFA+1,JXFE
      DFRACT = DBLE(IVALV(JXOP))
      NEXPU  = NEXPV(JXOP) + NEXPM
      IF (NEXPU.EQ.0)              GO TO 69
      IF (NEXPU.LT.0)              GO TO 68
      IF (NEXPU.GT.9)              GO TO 68
      DO 67  J=1,NEXPU
   67 DFRACT = DFRACT * 10.
      GO TO 69

   68 DFRACT = DFRACT * 10.D0**NEXPU
   69 DNUM   = DNUM + DFRACT
   70 IF (NEGM.NE.0)  DNUM = -DNUM
      GO TO 30

C--------          Reading octal or binary or hexadecimal

   71 J  = INDEX ('0OoBbXx', CHV(JJ+1))
      IF (J.EQ.0)                  GO TO 94
      JJ = JJ + 2
      IF (J.GE.6)                  GO TO 82
      IF (J.GE.4)                  GO TO 76

C--       octal

   72 IF (JJ.GT.JR)          GO TO 87
#if defined(CERNLIB_QASCII)
      K = ICHAR (CHV(JJ))
      IF (K.EQ.32)           GO TO 87
      K = K - 48
      IF (K.LT.0)            GO TO 86
      IF (K.GE.8)            GO TO 86
#endif
#if defined(CERNLIB_QEBCDIC)
      K = ICHAR (CHV(JJ))
      IF (K.EQ.64)           GO TO 87
      K = K - 240
      IF (K.LT.0)            GO TO 86
      IF (K.GE.8)            GO TO 86
#endif
#if !defined(CERNLIB_QISASTD)
      NUM(1) = ISHFTL(NUM(1),3) + K
#endif
#if defined(CERNLIB_QISASTD)
      NUM(1) = ISHFT(NUM(1),3) + K
#endif
      NDIGT = NDIGT + 1
      JJ    = JJ + 1
      GO TO 72

C--       binary

   76 IF (JJ.GT.JR)          GO TO 87
#if defined(CERNLIB_QASCII)
      K = ICHAR (CHV(JJ))
      IF (K.EQ.32)           GO TO 87
      K = K - 48
      IF (K.LT.0)            GO TO 86
      IF (K.GE.2)            GO TO 86
#endif
#if defined(CERNLIB_QEBCDIC)
      K = ICHAR (CHV(JJ))
      IF (K.EQ.64)           GO TO 87
      K = K - 240
      IF (K.LT.0)            GO TO 86
      IF (K.GE.2)            GO TO 86
#endif
#if !defined(CERNLIB_QISASTD)
      NUM(1) = ISHFTL(NUM(1),1) + K
#endif
#if defined(CERNLIB_QISASTD)
      NUM(1) = ISHFT(NUM(1),1) + K
#endif
      NDIGT = NDIGT + 1
      JJ    = JJ + 1
      GO TO 76

C--       hexadecimal

   82 IF (JJ.GT.JR)          GO TO 87
#if defined(CERNLIB_QASCII)
      K = ICHAR (CHV(JJ))
      IF (K.EQ.32)           GO TO 87
      K = K - 48
      IF (K.LT.0)            GO TO 86
      IF (K.LT.10)           GO TO 84
      K = K - 7
      IF (K.LT.8)            GO TO 86
      IF (K.LT.16)           GO TO 84
      K = K - 32
      IF (K.LT.8)            GO TO 86
      IF (K.GE.16)           GO TO 86
#endif
#if defined(CERNLIB_QEBCDIC)
      K = ICHAR (CHV(JJ))
      IF (K.EQ.64)           GO TO 87
      K = K - 240
      IF (K.GE.10)           GO TO 86
      IF (K.GE.0)            GO TO 84
      K = K + 57
      IF (K.GE.16)           GO TO 86
      IF (K.GE.8)            GO TO 84
      K = K + 64
      IF (K.GE.16)           GO TO 86
      IF (K.LT.8)            GO TO 86
#endif
#if !defined(CERNLIB_QISASTD)
   84 NUM(1) = ISHFTL(NUM(1),4) + K
#endif
#if defined(CERNLIB_QISASTD)
   84 NUM(1) = ISHFT(NUM(1),4) + K
#endif
      NDIGT = NDIGT + 1
      JJ    = JJ + 1
      GO TO 82

   86 JTERM = -1
   87 MODE  = 1
      IF (NDIGT.EQ.0)              GO TO 91
      GO TO 30

C----              Special error exits

   94 NERR = -1
   93 NERR = NERR - 1
   92 NERR = NERR - 1
   91 NERR = NERR - 1
   90 MODE = NERR
      NESLAT = JJ
      NGSLAT = JJ
      RETURN
      END
