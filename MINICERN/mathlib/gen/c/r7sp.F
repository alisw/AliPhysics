*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/04/01 15:01:57  mclareni
* Mathlib gen
*
*
#include "gen/pilot.h"
#if defined(CERNLIB_SINGLE)
      SUBROUTINE C309R7(A,B,IBEG,INUM,XX,EPS)
C
C*******************************************************************
C
C  RCF converts polynomial A to the corresponding continued
C         fraction, in 'normal'  form with coefficients B
C         by the 'P algorithm' of Patry & Gupta
C
C   A(z) = A1/z + A2/z**3 + A3/z**5 + ... + An/z**(2n-1)
C
C   B(z) = B1/z+ B2/z+ B3/z+ .../(z+ Bn/z)
C
C  data:
C   A     vector A(k), k=1,INUM         input
C   B     vector B(k), k=IBEG,INUM      output
C   IBEG  order of first coef. calc.    input
C   INUM  order of A, even or odd       input
C   XX    auxiliary vector of length .ge. length of vector B
C         caller provides space for A,B,XX
C   Note that neither of the first two terms A(1) A(2) should be zero
C          & the user can start the calculation with any value of
C          IBEG provided the c.f. coefs have been already
C          calculated up to INUM = IBEG-1
C          & the method breaks down as soon as the absolute value
C          of a c.f. coef. is less than EPS.    At the time of the
C          break up  INUM has been replaced by minus times the number
C          of this coefficient.
C   algorithm: J. Patry & S. Gupta, EIR-Bericht 247, November 1973
C              Eidg. Institut fur Reaktorforschung
C              Wuerenlingen, Switzerland
C   see also:  Haenggi, Roesel & Trautmann,
C              J. Comput. Phys., v. 137, (1980) 242-258
C   note:      restart procedure modified by I.J.Thompson
C
C*******************************************************************
C
      IMPLICIT COMPLEX(A-H,O-Z)
      DIMENSION A(100),B(100),XX(2,100)
      LOGICAL EVEN
      REAL EPS

      IBN=INUM
      IF(IBEG .GT. 4) GO TO 50
      IF(IBEG .EQ. 4) GO TO 20
      B(1)=A(1)
      IF(IBN .GE. 2) B(2)=-A(2)/A(1)
      IF(IBN .LT. 3) RETURN
      X0=A(3)/A(2)
      XX(2,1)=B(2)
      XX(1,1)=-X0
      XX(1,2)=0
      B(3)=-X0-B(2)
      X0=-B(3)*A(2)
      M=3
      MP12=2
      EVEN=.TRUE.
      IF(IBN .LE. 3) RETURN
   20 IF(ABS(B(3)) .LT. EPS*ABS(X0)) THEN
       INUM=-M
       RETURN
      END IF
      M=4
   30 X1=A(M)
      M2M1=MP12
      MP12=M2M1+1
      IF(EVEN) MP12=M2M1
      DO 40 K = 2,MP12
   40 X1=X1+A(M-K+1)*XX(1,K-1)
      B(M)=-X1/X0
      IF(M .GE. IBN) RETURN
   50 IF(ABS(B(M)) .LT. EPS*ABS(X0)) THEN
       INUM=-M
       RETURN
      END IF
      DO 60 K = M2M1,2,-1
   60 XX(2,K)=XX(1,K)+B(M)*XX(2,K-1)
      XX(2,1)=XX(1,1)+B(M)
      DO 70 K = 1,M2M1
      X0=XX(2,K)
      XX(2,K)=XX(1,K)
   70 XX(1,K)=X0
      X0=X1
      XX(1,M2M1+1)=0
      M=M+1
      EVEN=.NOT.EVEN
      GO TO 30
      END
      FUNCTION C309R8(Z,ACC)
      COMPLEX C309R8,Z
      REAL ACC
      REAL X,Y,AX,AY,A

      X=REAL(Z)
      Y=AIMAG(Z)
      AX=ABS(X)
      AY=ABS(Y)
      A=5*ACC*(AX+AY)
      IF(AX .LT. A) X=0
      IF(AY .LT. A) Y=0
      C309R8=CMPLX(X,Y)
      RETURN
      END
#endif
