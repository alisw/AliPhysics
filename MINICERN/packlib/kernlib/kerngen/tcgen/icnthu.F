*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/02/15 17:49:46  mclareni
* Kernlib
*
*
#include "kerngen/pilot.h"
      FUNCTION ICNTHU (TEXT,POSS,NPOSS)
C
C CERN PROGLIB# M432    ICNTHU          .VERSION KERNFOR  4.36  930602
C ORIG. 04/10/88, JZ
C
C-    Match TEXT against upper case POSS(NPOSS), case insensitive

      CHARACTER*(*)  TEXT,POSS(99)
      CHARACTER*1    CHP, CHT

#include "kerngen/qnatch.inc"
* Ignoring t=pass

      NTX = LEN(TEXT)
      NPO = LEN(POSS(1))
      JPOSS  = 0

   11 JPOSS = JPOSS + 1
      IF (JPOSS.GT.NPOSS)    GO TO 98
      JC = 1
      CHP = POSS(JPOSS)(1:1)
   12 CHT = TEXT(JC:JC)
      IF (CHT.NE.CHP)  THEN
          IVT = ICHAR(CHT)
          IVP = ICHAR(CHP)
#if defined(CERNLIB_QASCII)
          IF (IVT-IVP.NE.32) GO TO 11
          IF (IVT.LT.97)     GO TO 11
          IF (IVT.GT.122)    GO TO 11
#endif
#if defined(CERNLIB_QEBCDIC)
          IF (IVP-IVT.NE.64)      GO TO 11
          IF (NATCH(IVT+1).NE.3)  GO TO 11
#endif
        ENDIF
      JC = JC + 1
      IF (JC.GT.NPO)         GO TO 99
      CHP = POSS(JPOSS)(JC:JC)
      IF (CHP.EQ.'*')        GO TO 41
      IF (CHP.EQ.' ')        GO TO 31
      IF (JC.LE.NTX)         GO TO 12
      GO TO 11

C--                POSS blank terminated

   31 IF (JC.GT.NTX)         GO TO 99
      CHT = TEXT(JC:JC)
      IF (CHT.EQ.' ')        GO TO 99
      GO TO 11

C--                Check continuation after '*'

   41 IF (JC.GT.NTX)         GO TO 99
      CHT = TEXT(JC:JC)
      IF (CHT.EQ.' ')        GO TO 99
      JC = JC + 1
      IF (JC.GT.NPO)         GO TO 99
      CHP = POSS(JPOSS)(JC:JC)
      IF (CHP.EQ.' ')        GO TO 11
      IF (CHP.EQ.'*')        GO TO 99
      IF (CHT.EQ.CHP)        GO TO 41
      IVT = ICHAR(CHT)
      IVP = ICHAR(CHP)
#if defined(CERNLIB_QASCII)
      IF (IVT-IVP.NE.32)     GO TO 11
      IF (IVT.LT.97)         GO TO 11
      IF (IVT.GT.122)        GO TO 11
#endif
#if defined(CERNLIB_QEBCDIC)
      IF (IVP-IVT.NE.64)     GO TO 11
      IF (NATCH(IVT+1).NE.3) GO TO 11
#endif
      GO TO 41

   98 JPOSS = 0
   99 ICNTHU = JPOSS
      RETURN
      END
