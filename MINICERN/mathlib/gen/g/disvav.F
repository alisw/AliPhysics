*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/04/01 15:02:47  mclareni
* Mathlib gen
*
*
#include "gen/pilot.h"
      FUNCTION DISVAV(X,I)
C
C     DISVAV(X,I) COMPUTES FOR I=0 THE VALUE OF THE DENSITY, FOR I=1
C     THE VALUE OF THE CUMULATIVE DISTRIBUTION FUNCTION OF THE
C     VAVILOV DISTRIBUTION AT THE POINT X
C
      COMMON /VAVILI/ T0,T1,T,OMEGA
      COMMON /VAVILO/ A(155),B(155),N
C
      IF(X .LT. T0) GO TO 3
      IF(X .GT. T1) GO TO 4
C
      Y=X-T0
C
      TRR=0.
      IF(I .EQ. 1) TRR=Y/T
C
      Z=OMEGA*Y-3.1415926535898
      COF=2.0*COS(Z)
      A1=0.
      A0=A(1)
      NN=N-1
      DO 1 K = 2,N
      A2=A1
      A1=A0
    1 A0=A(K)+COF*A1-A2
C
      TRR=TRR+0.5*(A0-A2)
      A1=0.
      A0=B(1)
      DO 2 K = 2,NN
      A2=A1
      A1=A0
    2 A0=B(K)+COF*A1-A2
C
      DISVAV=TRR+A0*SIN(Z)
      RETURN
C
    3 DISVAV=0.0
      RETURN
C
    4 DISVAV=I
      RETURN
C
      END
