*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/02/15 17:49:54  mclareni
* Kernlib
*
*
#include "kerngen/pilot.h"
      SUBROUTINE TRCHUL (A,B,N)
C
C CERN PROGLIB# F112    TRCHUL          .VERSION KERNFOR  4.16  870601
C ORIG. 18/12/74 WH
C
#if defined(CERNLIB_INTDOUBL)
      DOUBLE PRECISION SUM, R, DC
#endif
      DIMENSION A(*),B(*)
C
      KPIV = (N*N+N)/2
C
      I = N
   10 IPIV = KPIV
      R = A(IPIV)
C
   20 SUM = 0.
      IF (I.EQ.N)             GO TO 40
      IF (R.EQ.0.)            GO TO 42
      ID = IPIV
      KD = KPIV
      NSTEP = I
C
   30 KD = KD + NSTEP
      ID = ID + NSTEP
      NSTEP = NSTEP + 1
      SUM = SUM + B(ID)*B(KD)
      IF (NSTEP.LT.N)         GO TO 30
C
   40 SUM = A(KPIV) - SUM
   42 IF (KPIV.LT.IPIV)       GO TO 50
#if !defined(CERNLIB_INTDOUBL)
      DC =  SQRT (SUM)
#endif
#if defined(CERNLIB_INTDOUBL)
      DC = DSQRT (SUM)
#endif
      B(KPIV) = DC
      IF (R .GT. 0.0) R = 1./DC
        GO  TO                60
C
   50 B(KPIV) = SUM*R
C
   60 KPIV = KPIV - 1
      IF (KPIV.GT.IPIV-I)     GO TO 20
C
      I = I - 1
      IF (I.GT.0)             GO TO 10
C
      RETURN
      END
