*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/02/15 17:49:53  mclareni
* Kernlib
*
*
#include "kerngen/pilot.h"
      SUBROUTINE TLRES (A,B,AUX)
C
C CERN PROGLIB# E230    TLRES           .VERSION KERNFOR  2.09  751101
C ORIG. 11/05/74 WH+WM
C
C.  SUBROUTINE TLRES          L.S. REDIDUAL VECTOR          HART/MATT
C.
C.        INVERSE HOUSEHOLDER TRANSFORMATION APPLIED TO LSQ RESIDUALS
C.        THE RESIDUALS ARE RETURNED IN B
C.  ARGUMENTS
C.        A,B,AUX,(M1,M,N,L,IER) WHICH ARE DEFINED AS FOR TLS
C.  REMARK
C.        CONSTRAINED AND UNCONSTRAINED VERSION COMBINED
C.
C.-------------------------------------------------------------------
C
      COMMON /TLSDIM/ M1,M,N,L,IER
      COMMON /SLATE/ BETA,I,IB,IST,KK,KN,K1,LN,LV,KKEND,DUM(30)
      DIMENSION A(*), B(*), AUX(*)
C
      IF(IER.EQ.0) GO TO 100
      IER = IABS(IER)
C
C     SET FIRST IER ELEMENTS OF EACH B VECTOR TO ZERO
      LN = L * IER
      DO 10 I=1,LN
   10 B(I)=0.
C
C      TRANSFORMATION LOOP
      K1 = MAX (N,L)
      IST = (IER-1) * (N+1) + 1
      KKEND = IER - M1
      DO 30 KK=1,KKEND
      LV = M - IER + KK
      KN = K1 + IER - KK + 1
      BETA = -1./(AUX(KN)*A(IST))
      IB = (IER-KK) * L + 1
      IF (LV.EQ.1)                     GO TO     20
      CALL TLSTEP(A(IST),B(IB),N,L,LV,L,BETA)
      GO TO        30
   20 DO           25        J=1,L
      JST = IB + J - 1
   25 B(JST) = B(JST)*(1.-BETA*A(IST)**2)
   30 IST = IST -N -1
  100 RETURN
      END
