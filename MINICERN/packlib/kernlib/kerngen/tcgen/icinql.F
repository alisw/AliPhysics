*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/02/15 17:49:45  mclareni
* Kernlib
*
*
#include "kerngen/pilot.h"
      FUNCTION ICINQL (TEXT,POSS,NPOSS)
C
C CERN PROGLIB# M432    ICINQL          .VERSION KERNFOR  4.36  930602
C ORIG. 24/05/93, JZ
C
C-    Match TEXT against lower case POSS(NPOSS), case insensitive

      CHARACTER*(*)  TEXT, POSS(99)
      CHARACTER*1    CHT,  CHP

#include "kerngen/qnatch.inc"
* Ignoring t=pass

      NTX = LEN(TEXT)
      NPO = LEN(POSS(1))

      JPOSS = 0
   11 JPOSS = JPOSS + 1
      IF (JPOSS.GT.NPOSS)    GO TO 98
      CHP = POSS(JPOSS)(1:1)
      CHT = TEXT(1:1)
      JC  = 1
   12 IF (CHT.NE.CHP)  THEN
          IVT = ICHAR(CHT)
          IVP = ICHAR(CHP)
#if defined(CERNLIB_QASCII)
          IF (IVP-IVT.NE.32) GO TO 11
          IF (IVT.LT.65)     GO TO 11
          IF (IVT.GT.90)     GO TO 11
#endif
#if defined(CERNLIB_QEBCDIC)
          IF (IVT-IVP.NE.64)      GO TO 11
          IF (NATCH(IVT+1).NE.4)  GO TO 11
#endif
        ENDIF
      JC = JC + 1
      IF (JC.GT.NPO)         GO TO 99
      CHP = POSS(JPOSS)(JC:JC)
      IF (CHP.EQ.' ')        GO TO 71
      IF (JC.GT.NTX)         GO TO 11
      CHT = TEXT(JC:JC)
      IF (CHT.NE.'*')        GO TO 12

C----         First '*' seen on TEXT

   61 JC = JC + 1
      IF (JC.GT.NTX)         GO TO 99
      CHT = TEXT(JC:JC)
      IF (CHT.EQ.'*')        GO TO 99
      IF (CHT.EQ.' ')        GO TO 99
      IF (CHP.NE.CHT)  THEN
          IVT = ICHAR(CHT)
          IVP = ICHAR(CHP)
#if defined(CERNLIB_QASCII)
          IF (IVP-IVT.NE.32) GO TO 11
          IF (IVT.LT.65)     GO TO 11
          IF (IVT.GT.90)     GO TO 11
#endif
#if defined(CERNLIB_QEBCDIC)
          IF (IVT-IVP.NE.64)     GO TO 11
          IF (NATCH(IVT+1).NE.4) GO TO 11
#endif
        ENDIF
      IF (JC.GT.NPO)         GO TO 99
      CHP = POSS(JPOSS)(JC:JC)
      IF (CHP.NE.' ')        GO TO 61
      GO TO 99

C----         POSS blank terminated

   71 IF (JC.GT.NTX)         GO TO 99
      CHT = TEXT(JC:JC)
      IF (CHT.EQ.' ')        GO TO 99
      IF (CHT.EQ.'*')        GO TO 99
      GO TO 11

   98 JPOSS = 0
   99 ICINQL = JPOSS
      RETURN
      END
