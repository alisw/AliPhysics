*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/04/01 15:02:43  mclareni
* Mathlib gen
*
*
#include "gen/pilot.h"
C     This will be GAMDIS,IF=DOUBLE and GAMDIS64,IF=-DOUBLE
      FUNCTION GAMDIS(X,A)

C     Calculates the gamma distribution function
C
C        G(x,a) = (1/gamma(a)) * int(0,x)[exp(-t) * t**(a-1)] dt.
C
C     Based on
C        W. Gautschi, ALGORITHM 542 Incomplete Gamma Functions,
C        ACM Trans. Math. Software 5 (1979) 482-489.

      CHARACTER NAME*(*)
      CHARACTER*80 ERRTXT
      PARAMETER (NAME = 'GAMDIS')
      DIMENSION C(14)

      PARAMETER (EPS = 1E-5, EPS1 = 5E-7)
      PARAMETER (ALH = -0.69314 72)
      PARAMETER (Z1 = 1, HALF = Z1/2, QUAR = Z1/4)
      PARAMETER (C1 = 3*Z1/2, KMAX = 300)

      DATA C
     1/ 0.5772157,-0.6558781,-0.0420026, 0.1665386,-0.0421977,
     2 -0.0096220, 0.0072189,-0.0011652,-0.0002152, 0.0001281,
     3 -0.0000201,-0.0000013, 0.0000011,-0.0000002/

      HST=0
      IF(X .EQ. 0) GO TO 99
      IF(X .LT. 0 .OR. A .LE. 0) THEN
       WRITE(ERRTXT,101) X,A
       CALL MTLPRT(NAME,'G106.1',ERRTXT)
       GO TO 99
      ELSE
       ALX=LOG(X)
      ENDIF
      IF(X .LT. QUAR) THEN
       ALFA=ALH/ALX
      ELSE
       ALFA=X+QUAR
      ENDIF
      IF(A .GT. ALFA) THEN
       TERM=1
       SUM=1
       DO 1 K = 1,KMAX
       TERM=X*TERM/(A+K)
       SUM=SUM+TERM
       IF(ABS(TERM) .LE. EPS*SUM) GO TO 2
    1  CONTINUE
       GO TO 98
    2  HST=SUM*EXP(A*ALX-X-ALOGAM(1+A))
      ELSEIF(X .GT. C1) THEN
       P=0
       S=1-A
       Q=(X+S)*(X-1-A)
       R=4*(X+S)
       TERM=1
       SUM=1
       RHO=0
       DO 3 K = 2,KMAX
       P=P+S
       Q=Q+R
       R=R+8
       S=S+2
       T=P*(1+RHO)
       RHO=T/(Q-T)
       TERM=RHO*TERM
       SUM=SUM+TERM
       IF(ABS(TERM) .LE. EPS*SUM) GO TO 4
    3  CONTINUE
       GO TO 98
    4  HST=1-(A*SUM/(X+1-A))*EXP(A*ALX-X-ALOGAM(1+A))
      ELSE
       IF(A .LT. HALF) THEN
        SUM=C(14)
        DO 12 K = 13,1,-1
   12   SUM=A*SUM+C(K)
        GA=-SUM/(1+A*SUM)
        Y=A*ALX
        IF(ABS(Y) .GE. 1) THEN
         U=GA-(EXP(Y)-1)/A
        ELSE
         SUM=1
         TERM=1
         DO 7 K = 2,KMAX
         TERM=Y*TERM/K
         SUM=SUM+TERM
         IF(ABS(TERM) .LE. EPS1*SUM) GO TO 8
    7    CONTINUE
         GO TO 98
    8    U=GA-SUM*ALX
        ENDIF
       ELSE
        U=GAMMA(A)-EXP(A*ALX)/A
       ENDIF
       P=A*X
       Q=A+1
       R=A+3
       TERM=1
       SUM=1
       DO 9 K = 2,KMAX
       P=P+X
       Q=Q+R
       R=R+2
       TERM=-P*TERM/Q
       SUM=SUM+TERM
       IF(ABS(TERM) .LE. EPS1*SUM) GO TO 10
    9  CONTINUE
       GO TO 98
   10  HST=1-A*(U+SUM*EXP((1+A)*ALX)/(1+A))/GAMMA(1+A)
      ENDIF
   99 GAMDIS=HST
      RETURN

   98 WRITE(ERRTXT,102) X,A
      CALL MTLPRT(NAME,'G106.2',ERRTXT)
      GO TO 99
  101 FORMAT('ILLEGAL ARGUMENT(S) X = ',E15.8,'  A = ',E15.8)
  102 FORMAT('PROBLEMS WITH CONVERGENCE, X = ',E15.8,'  A = ',E15.8)
      END
