*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/02/15 17:49:54  mclareni
* Kernlib
*
*
#include "kerngen/pilot.h"
      SUBROUTINE TRCHLU (A,B,N)
C
C CERN PROGLIB# F112    TRCHLU          .VERSION KERNFOR  4.16  870601
C ORIG. 18/12/74 W.HART
C
#if defined(CERNLIB_INTDOUBL)
      DOUBLE PRECISION SUM, R, DC
#endif
      DIMENSION A(*),B(*)
C
      IPIV = 0
C
      I = 0
   10 I = I + 1
      IPIV = IPIV + I
      KPIV = IPIV
      R = A(IPIV)
C
      DO 60 J=I,N
      SUM = 0.
      IF (I.EQ.1)             GO TO 40
      IF (R.EQ.0.)            GO TO 42
      ID = IPIV - I + 1
      KD = KPIV - I + 1
C
   30 SUM = SUM + B(KD)*B(ID)
      KD = KD + 1
      ID = ID + 1
      IF (ID.LT.IPIV)         GO TO 30
C
   40 SUM = A(KPIV) - SUM
   42 IF (J.NE.I)             GO TO 50
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
   60 KPIV = KPIV + J
C
      IF (I.LT.N)             GO TO 10
C
      RETURN
      END
