*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/04/01 15:02:17  mclareni
* Mathlib gen
*
*
#include "gen/pilot.h"
#if !defined(CERNLIB_DOUBLE)
      SUBROUTINE RDEQMR(N,XA,XZ,Y,H0,EPS,SUB,W)
#endif
#if defined(CERNLIB_DOUBLE)
      SUBROUTINE DDEQMR(N,XA,XZ,Y,H0,EPS,SUB,W)
#include "gen/imp64.inc"
#endif
C     Based on a modification of the Runge-Kutta method suggested
C     by Merson. See G.N. Lance, Numerical Methods for High speed
C     Computers, Iliffe & Sons, London 1960, pp. 56-57

      CHARACTER NAME*(*)
      CHARACTER*80 ERRTXT
#if !defined(CERNLIB_DOUBLE)
      PARAMETER (NAME = 'RDEQMR')
#endif
#if defined(CERNLIB_DOUBLE)
      PARAMETER (NAME = 'DDEQMR')
#endif
      LOGICAL LER,LFN

      DIMENSION Y(*),W(N,*)

      PARAMETER (DELTA = 1D-14)
      PARAMETER (Z1 = 1, R2 = Z1/2, R3 = Z1/3)
      PARAMETER (R4 = 3*Z1/8, R5 = 3*Z1/2, R6 = 9*Z1/2)
      PARAMETER (R7 = 4*Z1/3, R0 = Z1/32)

#if !defined(CERNLIB_DOUBLE)
      ENTRY DEQMR(N,XA,XZ,Y,H0,EPS,SUB,W)
#endif

      IF(N .LT. 1 .OR. XA .EQ. XZ .OR. H0 .EQ. 0) RETURN
      DELTAX=DELTA*ABS(XZ-XA)
      EPS5=5*ABS(EPS)
      EPS0=R0*EPS5
      X=XA
      H1=SIGN(ABS(H0),XZ-XA)
      SGH=SIGN(Z1,H1)

   12 IF(SGH*(X+H1-XZ) .LT. 0) THEN
       HH=H1
       H0=H1
       LFN=.FALSE.
      ELSE
       HH=XZ-X
       IF(ABS(HH) .LT. DELTAX) THEN
        DO 10 I = 1,N
   10   Y(I)=W(I,6)
        RETURN
       END IF
       LFN=.TRUE.
      END IF
      S2=R2*HH
      S3=R3*HH
      S7=R7*HH
      X1=X+HH
      X2=X+S2
      X3=X+S3

      CALL SUB(X,Y,W(1,1))
      DO 1 I = 1,N
      W(I,1)=S3*W(I,1)
    1 W(I,6)=Y(I)+W(I,1)

      CALL SUB(X3,W(1,6),W(1,2))
      DO 2 I = 1,N
      W(I,2)=S3*W(I,2)
    2 W(I,6)=Y(I)+R2*(W(I,1)+W(I,2))

      CALL SUB(X3,W(1,6),W(1,3))
      DO 3 I = 1,N
      W(I,3)=S3*W(I,3)
      W(I,2)=3*W(I,3)
    3 W(I,6)=Y(I)+R4*(W(I,1)+W(I,2))

      CALL SUB(X2,W(1,6),W(1,4))
      DO 4 I = 1,N
      W(I,4)=S7*W(I,4)
    4 W(I,6)=Y(I)+R5*(W(I,1)-W(I,2)+W(I,4))

      CALL SUB(X1,W(1,6),W(1,5))
      DO 5 I = 1,N
      W(I,5)=S3*W(I,5)
    5 W(I,6)=Y(I)+R2*(W(I,1)+W(I,4)+W(I,5))

      DO 8 I = 1,N
      W(I,2)=ABS(W(I,1)-R6*W(I,3)+W(I,4)-R2*W(I,5))
      W(I,1)=ABS(W(I,6))
      IF(W(I,2) .GT. EPS5*W(I,1)) THEN
       H1=R2*HH
       IF(ABS(H1) .LT. DELTAX) THEN
        WRITE(ERRTXT,101) X
        CALL MTLPRT(NAME,'D202.1',ERRTXT)
        RETURN
       END IF
       GO TO 12
      END IF
    8 CONTINUE
      LER=.TRUE.
      DO 7 I = 1,N
    7 LER=LER .AND. W(I,2) .LT. EPS0*W(I,1)
      DO 9 I = 1,N
    9 Y(I)=W(I,6)
      IF(LER) THEN
       H0=H1+H1
       H1=HH+HH
      END IF
      IF(LFN) RETURN
      X=X1
      GO TO 12
  101 FORMAT('TOO HIGH ACCURACY REQUIRED NEAR  X = ',1P,D15.8)
      END
