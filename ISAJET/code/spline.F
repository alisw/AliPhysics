#include "isajet/pilot.h"
      SUBROUTINE SPLINE(X,C,N,IBCBEG,IBCEND)
C**********************************************************************
C* Computes the coefficient of a cubic interpolating spline. The X(i),*
C* i=1,...,N, are the knots or x-values of data points; C(1,i) are    *
C* the corresponding y-values. N is the number of data points (>3!).  *
C* IBCBEG = 0 means that the slope at X(1) is unknown, in which case  *
C* it is determined from requiring a smooth 3rd derivative at x(2);   *
C* IBCBEG = 1 means that the slope is known, in which case it has to  *
C* be stored in C(1,2). IBCEND has the same meaning for the end of    *
C* x-region; if IBCEND = 1, the slope is to be stored in C(2,N). The  *
C* routine then computes the coefficients C(l,i) of the i-th spline,  *
C* written in the form                                                *
C* f_i(x) = C(1,i) + h_i C(2,i) + h_i^2 C(3,i) + h_i^3 C(4,i), where  *
C* h_i = x - X(I).                                                    *
C  Modified from contributed subroutine by M. Drees, 1/14/99
C**********************************************************************
#if defined(CERNLIB_IMPNONE)
      IMPLICIT NONE
#endif
      INTEGER N,IBCBEG,IBCEND,I,L,M,J
      REAL C(4,N),X(N),G,DTAU,DIVDF1,DIVDF3
C
      L = N - 1
      DO 10 M = 2, N
      C(3,M) = X(M) - X(M-1)
   10 C(4,M) = (C(1,M) - C(1,M-1))/C(3,M)
C          First slope unknown
      IF(IBCBEG.EQ.0) THEN
        C(4,1) = C(3,3)
        C(3,1) = C(3,2) + C(3,3)
        C(2,1) = ( (C(3,2)+2.0*C(3,1))*C(4,2)*C(3,3) +
     &             C(3,2)**2*C(4,3) )/C(3,1)
C          First slope already known
      ELSE
        C(4,1) = 1.0
        C(3,1) = 0.0
      ENDIF
C          Forward pass of Gauss elimination
      DO 20 M = 2, L
        G = -C(3,M+1)/C(4,M-1)
        C(2,M) = G*C(2,M-1) + 3.0*(C(3,M)*C(4,M+1)+C(3,M+1)*C(4,M))
   20   C(4,M) = G*C(3,M-1) + 2.0*(C(3,M)+C(3,M+1))
 
      IF(IBCEND.EQ.0) THEN
        G = C(3,N-1) + C(3,N)
        C(2,N) = ( (C(3,N)+2.0*G)*C(4,N)*C(3,N-1) +
     &             C(3,N)**2*(C(1,N-1)-C(1,N-2))/C(3,N-1) )/G
        G = -G/C(4,N-1)
        C(4,N) = (G+1.0)*C(3,N-1) + C(4,N)
        C(2,N) = ( G*C(2,N-1) + C(2,N) )/C(4,N)
      ENDIF
C          Back substitution
      DO 30 J = L,1,-1
   30 C(2,J) = ( C(2,J) - C(3,J)*C(2,J+1) )/C(4,J)
C          Computation of coefficients
      DO 40 I = 2,N
        DTAU = C(3,I)
        DIVDF1 = (C(1,I)-C(1,I-1))/DTAU
        DIVDF3 = C(2,I-1) + C(2,I) - 2.0*DIVDF1
        C(3,I-1) = ( DIVDF1 - C(2,I-1) - DIVDF3 ) / DTAU
   40   C(4,I-1) = DIVDF3/DTAU/DTAU
C 
      RETURN
      END
