*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/04/01 15:02:25  mclareni
* Mathlib gen
*
*
#include "gen/pilot.h"
      FUNCTION TRISUM(X,A,N,B,M,IOP)
C
C     PURPOSE = TO COMPUTE THE SUM OF THE TRIGONOMETRIC SERIES
C
C     TRISUM  = A(1)/2+A(2)*COS(X)+....................A(N)*COS((N-1)*X)
C                     +B(1)*SIN(X)+....................B(M)*SIN(M*X)
C
C     PARAMETERS
C
C     X       = FUNCTION ARGUMENT
C     A       = ARRAY DIMENSIONED A(N),KEEPING THE COEFFICIENTS FOR THE
C               COSINE TERMS
C     N       = THE NUMBER OF COEFFICIENTS IN THE COSINE SERIE
C     B       = ARRAY DIMENSIONED B(M),KEEPING THE COEFFICIENTS FOR THE
C               SINE TERMS
C     M       = THE NUMBER OF COEFFICIENTS IN THE SINE SERIE
C     IOP     = OPTION NUMBER , IOP = 1 , THE GENERAL CASE
C                               IOP = 2 , THE EVEN CASE
C                               IOP = 3 , THE ODD CASE
C
      DIMENSION A(N) , B(M)
C
C*UL 1000 COF=2.*COS(X)
      COF=2.*COS(X)
      TRISUM=0.
C
      IF (IOP.EQ.3) GO TO 1030
C
C     COMPUTE SUM OF COSINE TERMS
C
      A2 =0.
      A1 =A2
      A0 =A(N)
C
      IF (N.EQ.1) GO TO 1020
C
      NN=N-1
      DO 1010 K=1,NN
      A2=A1
      A1=A0
      INDEX=NN+1-K
      A0=A(INDEX)+COF*A1-A2
 1010 CONTINUE
C
 1020 TRISUM=.5*(A0-A2)
C
 1030 IF (IOP.EQ.2) RETURN
C
C     COMPUTE SUM OF SINUS TERMS
C
      A2 =0.
      A1 =A2
      A0 =B(M)
C
      IF (M.EQ.1) GO TO 1050
C
      NN=M-1
      DO 1040 K=1,NN
      A2=A1
      A1=A0
      INDEX=NN+1-K
      A0=B(INDEX)+COF*A1-A2
 1040 CONTINUE
C
 1050 TRISUM=TRISUM+A0*SIN(X)
C
      RETURN
C
      END
