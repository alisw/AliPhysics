*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/02/15 17:49:54  mclareni
* Kernlib
*
*
#include "kerngen/pilot.h"
      SUBROUTINE TRLTA (U,A,B,M,N)
C
C CERN PROGLIB# F112    TRLTA           .VERSION KERNFOR  4.15  861204
C ORIG. 18/12/74 WH
C
#if defined(CERNLIB_INTDOUBL)
      DOUBLE PRECISION SUM
#endif
      DIMENSION A(*),U(*),B(*)
C
      IPIV = 0
      MX = M*N
      MXPN = MX + N
      IB = 0
C
      I = 0
   10 I = I + 1
      IPIV = IPIV + I
C
   20 IU = IPIV
      NSTEP = I
      IB = IB + 1
      IA = IB
C
      SUM = 0.
   30 SUM = SUM + A(IA)*U(IU)
      IA = IA + N
      IU = IU + NSTEP
      NSTEP = NSTEP + 1
      IF (IA.LE.MX)           GO TO 30
C
      B(IB) = SUM
      IF (IA.LT.MXPN)         GO TO 20
C
      IF (I.LT.M)             GO TO 10
C
      RETURN
      END
