*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/02/15 17:49:43  mclareni
* Kernlib
*
*
#include "kerngen/pilot.h"
      SUBROUTINE CSETVI (NINTP,INTV,NBIASP,CHV,JLP,JRP,NSTEPP,IFLSQP)
C
C CERN PROGLIB# M432    CSETVI          .VERSION KERNFOR  4.36  930602
C ORIG. 01/04/92, JZ
C
C-    Set NINT decimal integers INTV(J)+NBIAS into CHV(JL:JR)
C-    every NSTEP columns, each with one leading blank
C-    and right-justified within its field;
C-    squeeze multiple blanks if IFLSQ non-zero.

      COMMON /SLATE/ NDSLAT,NESLAT,NFSLAT,NGSLAT, DUMMY(36)

      DIMENSION    NINTP(9), NBIASP(9)
      DIMENSION    JLP(9), JRP(9),  NSTEPP(9), IFLSQP(9)
      DIMENSION    INTV(*)
      CHARACTER    CHV*(*)


      NINT  = NINTP(1)
      NBIAS = NBIASP(1)
      JL    = JLP(1)
      NSTEP = NSTEPP(1)
      LIMIT = JRP(1) - NSTEP + 1
      IFLSQ = IFLSQP(1)
      JWK   = JL - 1

      DO 24  JIN=1,NINT
      IF (JWK.GE.LIMIT)  THEN
          IF (IFLSQ.EQ.0)          GO TO 91
          IF (JIN.EQ.1)            GO TO 91
          CALL CSQMBL (CHV,JL,JWK)
          JWK = NESLAT - 1
          IF (JWK.GE.LIMIT)        GO TO 91
        ENDIF
      JVAL = INTV(JIN) + NBIAS
      CALL CSETDI (JVAL,CHV,JWK+2,JWK+NSTEP)
   24 JWK = JWK + NSTEP

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
