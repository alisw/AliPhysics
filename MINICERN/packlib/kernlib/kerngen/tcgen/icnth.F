*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/02/15 17:49:46  mclareni
* Kernlib
*
*
#include "kerngen/pilot.h"
      FUNCTION ICNTH (TEXT,POSS,NPOSS)
C
C CERN PROGLIB# M432    ICNTH           .VERSION KERNFOR  4.36  930602
C ORIG. 04/10/88, JZ
C
C-    Match TEXT against POSS(NPOSS), case sensitive

      CHARACTER*(*)  TEXT,POSS(99)
      CHARACTER*1    CHP, CHT


      NTX = LEN(TEXT)
      NPO = LEN(POSS(1))
      JPOSS  = 0
C     PRINT 9009, NTX,NPO, TEXT
C9009 FORMAT ('NTX,NPO=',2I4,' find >',A,'<')

   11 JPOSS = JPOSS + 1
      IF (JPOSS.GT.NPOSS)    GO TO 98
      JC = 1
C     PRINT 9011, JPOSS,POSS(JPOSS)
C9011 FORMAT ('JPOSS=',I3,' take >',A,'<')
      CHP = POSS(JPOSS)(1:1)
   12 CHT = TEXT(JC:JC)
C     PRINT 9012, CHT,CHP
C9012 FORMAT ('@12 CHT/CHP= ',A,1X,A)
      IF (CHT.NE.CHP)        GO TO 11
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
C     PRINT 9031, JC,CHT
C9031 FORMAT('at col. JC',I3,' >',A,'<')
      IF (CHT.EQ.' ')        GO TO 99
      GO TO 11

C--                Check continuation after '*'

   41 IF (JC.GT.NTX)         GO TO 99
      CHT = TEXT(JC:JC)
      IF (CHT.EQ.' ')        GO TO 99
      JC = JC + 1
      IF (JC.GT.NPO)         GO TO 99
      CHP = POSS(JPOSS)(JC:JC)
C     PRINT 9044, CHT,CHP
C9044 FORMAT ('@44 CHT/CHP= ',A,1X,A)
      IF (CHP.EQ.' ')        GO TO 11
      IF (CHP.EQ.'*')        GO TO 99
      IF (CHT.EQ.CHP)        GO TO 41
      GO TO 11

   98 JPOSS = 0
   99 ICNTH = JPOSS
      RETURN
      END
