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
      FUNCTION C309R3(AA,BB,Z,EPS,JMAX,RE,FPMAX,N,X)
C
C     evaluate the HYPERGEOMETRIC FUNCTION
C                                             i
C            F (AA,BB;;Z) = SUM  (AA)  (BB)  Z / i!
C           2 0              i       i     i
C
C     to accuracy EPS with at most JMAX terms.
C
C     if the terms start diverging,
C     the corresponding continued fraction is found by RCF
C     & evaluated progressively by Steed's method to obtain convergence.
C
      IMPLICIT COMPLEX(A-H,O-Z)
      DIMENSION X(JMAX,4)
      LOGICAL FINITE
      REAL EP,EPS,AT,ATL,ABSC,RE,FPMAX

      ABSC(W)=ABS(REAL(W))+ABS(AIMAG(W))
      NINTC(W)=NINT(REAL(W))
C
      RE=0
      X(1,1)=1
      SUM=X(1,1)
      ATL=ABSC(X(1,1))
      F=SUM
      D=1
      DF=SUM
      J=0
      EP=EPS*(10*JMAX)
      MA=-NINTC(AA)
      MB=-NINTC(BB)
      FINITE=ABS(ABS(REAL(AA))-MA) .LT. EP .AND. ABS(AIMAG(AA)) .LT. EP
     1  .OR. ABS(ABS(REAL(BB))-MB) .LT. EP .AND. ABS(AIMAG(BB)) .LT. EP
      IMAX=JMAX
      IF(FINITE .AND. MA .GE. 0) IMAX=MIN(MA+1,IMAX)
      IF(FINITE .AND. MB .GE. 0) IMAX=MIN(MB+1,IMAX)
      DO 10 I = 2,IMAX
      X(I,1)=X(I-1,1)*Z*(AA+I-2)*(BB+I-2)/(I-1)
      IF(ABSC(X(I,1)) .GT. FPMAX) THEN
       N=0
       C309R3=SUM
       IF(.NOT.FINITE) RE=AT/ABSC(SUM)
       RETURN
      END IF
      AT=ABSC(X(I,1))
      IF(J .EQ. 0) THEN
       SUM=SUM+X(I,1)
       IF(AT .LT. ABSC(SUM)*EPS) THEN
        N=I
        C309R3=SUM
        IF(.NOT.FINITE) RE=AT/ABSC(SUM)
        RETURN
       END IF
      END IF
      IF(FINITE) GO TO 10
      IF(J .GT. 0 .OR. AT .GT. ATL .OR. I .GE. JMAX-2) J=J+1
      IF(J .EQ. 0) GO TO 10
      CALL C309R7(X(1,1),X(1,2),J,I,X(1,3),EPS)
      IF(I .LT. 0) THEN
       N=0
       C309R3=SUM
       IF(.NOT.FINITE) RE=AT/ABSC(SUM)
       RETURN
      END IF
      DO 50 K = MAX(J,2),I
      D=1/(D*X(K,2)+1)
      DF=DF*D-DF
      F=F+DF
      IF(ABSC(DF) .LT. ABSC(F)*EPS .OR.
     1   DF .EQ. 0 .AND. F .EQ. 0 .AND. I .GE. 4) THEN
       N=K
       C309R3=F
       RE=ABSC(DF)/ABSC(F)
       RETURN
      END IF
   50 CONTINUE
      J=I
   10 ATL=AT
      IF(.NOT.FINITE) I=-JMAX
      N=I
      C309R3=SUM
      IF(.NOT.FINITE) RE=AT/ABSC(SUM)
      RETURN
      END
#endif
