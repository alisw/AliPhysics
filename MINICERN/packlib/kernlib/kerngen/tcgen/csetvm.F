*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/02/15 17:49:43  mclareni
* Kernlib
*
*
#include "kerngen/pilot.h"
      SUBROUTINE CSETVM (NTIMEP,NINCP,JVALP,CHV,JLP,JRP,NSTEPP,IFLSQP)
C
C CERN PROGLIB# M432    CSETVM          .VERSION KERNFOR  4.36  930602
C ORIG. 01/04/92, JZ
C
C-    Set NTIME decimal integers JVAL + n*NINC into CHV(JL:JR)
C-    every NSTEP columns, each with one leading blank
C-    and right-justified within its field;
C-    squeeze multiple blanks if IFLSQ non-zero.

      COMMON /SLATE/ NDSLAT,NESLAT,NFSLAT,NGSLAT, DUMMY(36)

      DIMENSION    NTIMEP(9), NINCP(9), JVALP(9)
      DIMENSION    JLP(9), JRP(9),  NSTEPP(9), IFLSQP(9)
      CHARACTER    CHV*(*)


      NTIME = NTIMEP(1)
      NINC  = NINCP(1)
      JVAL  = JVALP(1)
      JL    = JLP(1)
      NSTEP = NSTEPP(1)
      LIMIT = JRP(1) - NSTEP + 1
      IFLSQ = IFLSQP(1)
      JWK   = JL - 1

      DO 24  JIN=1,NTIME
      IF (JWK.GE.LIMIT)  THEN
          IF (IFLSQ.EQ.0)          GO TO 91
          IF (JIN.EQ.1)            GO TO 91
          CALL CSQMBL (CHV,JL,JWK)
          JWK = NESLAT - 1
          IF (JWK.GE.LIMIT)        GO TO 91
        ENDIF
      CALL CSETDI (JVAL,CHV,JWK+2,JWK+NSTEP)
      JVAL = JVAL + NINC
   24 JWK  = JWK + NSTEP

      IF (IFLSQ.NE.0)  THEN
          CALL CSQMBL (CHV,JL,JWK)
          JWK = NESLAT - 1
        ENDIF

      NESLAT = JWK
      NGSLAT = 0
      RETURN

C----              Errors

   91 NESLAT = JWK
      NGSLAT = JIN
      RETURN
      END
