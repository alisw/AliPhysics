*
* $Id$
*
* $Log$
* Revision 1.2  1996/04/18 16:13:05  mclareni
* Incorporate changes from J.Zoll for version 3.77
*
* Revision 1.1.1.1  1996/03/06 10:47:22  mclareni
* Zebra
*
*
#include "zebra/pilot.h"
      SUBROUTINE MZIOCR (IOW)

C-    Crack I/O characteristic ready for use, system called

#include "zebra/quest.inc"
#include "zebra/mzioc.inc"
C--------------    END CDE                             --------------
      EQUIVALENCE (JIO,IQUEST(1))
      DIMENSION    IOW(9)
      DIMENSION    NBITVA(4), NBITVB(4), NBITVC(7)
#if defined(CERNLIB_QMVDS)
      SAVE         NBITVA, NBITVB, NBITVC
#endif
      DATA  NBITVA / 32,    16,   10,   8 /
      DATA  NBITVB / 29,    14,    9,   7 /
      DATA  NBITVC / 26,    11,    6,   4, 2, 1, 1 /

#include "zebra/q_jbit.inc"
#include "zebra/q_jbyt.inc"


      NWFODN = 0
      JFOCUR = 0

      JTYPR  = IOW(1)
      IOW1   = JBYT (JTYPR,17,16)

C--                IMMEDIATE CHARACTERISTIC

      IF (IOW1.NE.0)               GO TO 21
      IF (JTYPR.GE.8)              GO TO 21
      MFO(1) = JTYPR
      MFO(2) = -1
      JFOEND = 2
      JFOREP = 2
      RETURN

C----              GENERAL

   21 JFOEND = 0
      JFOREP = 0
      JIO    = 1
      JTYPR  = JBYT (IOW1,1,3)
      JFLAG  = JBIT (IOW1,4)
      JCLASS = JBYT (IOW1,14,3)
      JFL12  = 0

      GO TO ( 101, 201, 301, 401, 501, 601, 991), JCLASS

C--                CLASS 0 : 'CT -T'

      JFL12 = JFLAG + 1
      JTYP  = JBYT (IOW1,5,3)
      IF (JTYP.NE.0)  THEN
          MFO(1) = JTYP
          MFO(2) = JBYT (IOW1,8,6)
          JFOEND = 2
        ENDIF

   24 IF (JTYPR.EQ.7)      GO TO 28
      MFO(JFOEND+1) = JTYPR
      MFO(JFOEND+2) = JFL12 - 2
      JFOEND = JFOEND + 2
      JFOREP = JFOEND
      RETURN

   28 JFOREP = JFOEND
      MFO(JFOEND+1) = 7
      MFO(JFOEND+2) = 0
      JFOEND = JFOEND + 2
      RETURN

C--                CLASS 1 OR 2 : 'CT ... CT -T'  OR  'CT ... CT *CT'

  101 CONTINUE
  201 JFL12 = JCLASS
      IF (JTYPR.NE.0)              GO TO 821
      JTYPR = JBYT (IOW1,5,3)
      JBT   = 8
      GO TO 831

C--                CLASS 3 : 'CT / *T'  OR  '/ CT *T'

  301 JTYP = JBYT (IOW1,5,3)
      IF (JTYP.NE.0)  THEN
          MFO(1) = JTYP
          MFO(2) = JBYT (IOW1,8,6)
          JFOEND = 2
          IF (JFLAG.EQ.0)  JFOREP = 2
        ENDIF

      MFO(JFOEND+1) = JTYPR
      MFO(JFOEND+2) = 0
      JFOEND = JFOEND + 2
      RETURN

C--                CLASS 4 :  'CT / CT CT CT' OR 'CT CT / CT CT'

  401 JFOREP = 2*(JFLAG+1)
      JFLAG  = 0

C--                CLASS 5 : '/ CT ... CT'

  501 IF (JTYPR.EQ.0)              GO TO 830
      MFO(1) = JTYPR
      JFOEND = 2
      GO TO 821

C--                CLASS 6 : 'CT ... CT / CT ... CT'

  601 JFOREP = 2*JBYT(IOW1,1,4)
      JFLAG  = 1

C----              COMMON UNPACKING FOR CLASSES 1, 2, 4, 5, 6

  821 JIO = 2
      DO  822  JBT=5,11,3
      JTYP = JBYT (IOW1,JBT,3)
      IF (JTYP.EQ.0)               GO TO 823
      MFO(JFOEND+1) = JTYP
  822 JFOEND = JFOEND + 2
  823 NGRU   = JFOEND/2

C--                UNPACK I/O WORD 2

      IF (JFLAG.EQ.0)  THEN
          NBT = NBITVA(NGRU)
        ELSE
          NBT = NBITVB(NGRU)
        ENDIF

      JFOEND = 0
      JBT    = 1
      IOWN   = IOW(2)
      DO  824  JL=1,NGRU
      MFO(JFOEND+2) = JBYT(IOWN,JBT,NBT)
      JFOEND = JFOEND + 2
  824 JBT    = JBT + NBT
      IF (JFLAG.EQ.0)              GO TO 839

C--                UNPACK I/O WORDS 3, 4, ...

  825 NGRU = JBYT(IOWN,30,3)
      IF (NGRU.EQ.0)               GO TO 839
      JIO  = JIO + 1
      IF (JIO.EQ.17)               GO TO 991
      IOWN = IOW(JIO)
      JBTT = 1
      JBTC = 3*NGRU + 1
      NBT  = NBITVC(NGRU)

      DO  826  JL=1,NGRU
      MFO(JFOEND+1) = JBYT (IOWN,JBTT,3)
      MFO(JFOEND+2) = JBYT (IOWN,JBTC,NBT)
      JBTT   = JBTT + 3
      JBTC   = JBTC + NBT
  826 JFOEND = JFOEND + 2
      GO TO 825

C----              SUB-CLASSES ZERO

  830 JBT = 5
  831 DO  834  JL=JBT,11,3
      JTYP = JBYT (IOW1,JL,3)
      IF (JTYP.EQ.0)               GO TO 839
      MFO(JFOEND+1) = JTYP
      MFO(JFOEND+2) = 0
  834 JFOEND = JFOEND + 2

  839 IF (JFL12.NE.0)              GO TO 24
      RETURN

C-------           TROUBLE

  991 IQUEST(1) = -1
      MFO(1) = 0
      MFO(2) = -1
      JFOEND = 2
      RETURN
      END
*      ==================================================
#include "zebra/qcardl.inc"
