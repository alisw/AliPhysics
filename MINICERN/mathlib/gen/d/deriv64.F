*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/04/01 15:02:19  mclareni
* Mathlib gen
*
*
#include "gen/pilot.h"
#if !defined(CERNLIB_DOUBLE)
      SUBROUTINE RDERIV(F,X,DELTA,DFDX,RERR)
#endif
#if defined(CERNLIB_DOUBLE)
      SUBROUTINE DDERIV(F,X,DELTA,DFDX,RERR)
#include "gen/imp64.inc"
#endif
C     Computes the derivative f'(x) of f(x) at x = X. Based on
C     H. Rutishauser, Ausdehnung des Rombergschen Prinzips
C     (Extension of Romberg's Principle), Numer. Math. 5 (1963) 48-54


      DIMENSION DX(0:9),W(0:9,3),T(0:9,0:9),A(0:9)
      LOGICAL LEV(0:9),LMT

      CHARACTER NAME*(*)
      CHARACTER*80 ERRTXT
#if !defined(CERNLIB_DOUBLE)
      PARAMETER (NAME = 'RDERIV')
      PARAMETER (EPS = 5D-12)
#endif
#if defined(CERNLIB_DOUBLE)
      PARAMETER (NAME = 'DDERIV')
      PARAMETER (EPS = 5D-14)
#endif
      PARAMETER (Z1 = 1, S = Z1/10)

      DATA DX /0.0256D0, 0.0192D0, 0.0128D0, 0.0096D0, 0.0064D0,
     1         0.0048D0, 0.0032D0, 0.0024D0, 0.0016D0, 0.0012D0/

      DATA (LEV(K),K=0,8,2) /5*.TRUE./
      DATA (LEV(K),K=1,9,2) /5*.FALSE./

      DATA W(1,1) /1.33333 33333 333333D+00/
      DATA W(3,1) /1.06666 66666 666667D+00/
      DATA W(5,1) /1.01587 30158 730159D+00/
      DATA W(7,1) /1.00392 15686 274510D+00/

      DATA W(2,1) /3.33333 33333 333333D-01/
      DATA W(4,1) /6.66666 66666 666667D-02/
      DATA W(6,1) /1.58730 15873 015873D-02/
      DATA W(8,1) /3.92156 86274 509804D-03/

      DATA W(0,2) /2.28571 42857 142857D+00/
      DATA W(2,2) /1.16363 63636 363636D+00/
      DATA W(4,2) /1.03643 72469 635628D+00/
      DATA W(6,2) /1.00886 69950 738916D+00/
      DATA W(8,2) /1.00220 21042 329337D+00/

      DATA W(1,2) /1.28571 42857 142857D+00/
      DATA W(3,2) /1.63636 36363 636364D-01/
      DATA W(5,2) /3.64372 46963 562753D-02/
      DATA W(7,2) /8.86699 50738 916256D-03/
      DATA W(9,2) /2.20210 42329 336922D-03/

      DATA W(0,3) /1.80000 00000 000000D+00/
      DATA W(2,3) /1.12500 00000 000000D+00/
      DATA W(4,3) /1.02857 14285 714286D+00/
      DATA W(6,3) /1.00699 30069 930070D+00/
      DATA W(8,3) /1.00173 91304 347826D+00/

      DATA W(1,3) /8.00000 00000 000000D-01/
      DATA W(3,3) /1.25000 00000 000000D-01/
      DATA W(5,3) /2.85714 28571 428571D-02/
      DATA W(7,3) /6.99300 69930 069930D-03/
      DATA W(9,3) /1.73913 04347 826087D-03/

#if !defined(CERNLIB_DOUBLE)
      ENTRY DERIV(F,X,DELTA,DFDX,RERR)
#endif

      DEL=10*ABS(DELTA)
      IS=10

    4 IS=IS-1
      DEL=S*DEL
      IF(IS .EQ. 0 .OR. X+DEL*DX(9) .EQ. X) THEN
       DELTA=DEL
       DFDX=0
       RERR=1
       WRITE(ERRTXT,101) X
       CALL MTLPRT(NAME,'D401.1',ERRTXT)
       RETURN
      ENDIF
      DO 1 K = 0,9
      H=DEL*DX(K)
      T(K,0)=(F(X+H)-F(X-H))/(H+H)
    1 A(K)=T(K,0)

      IF(A(0) .GE. A(9)) THEN
       DO 5 K = 0,9
    5  A(K)=-A(K)
      ENDIF

      LMT=.TRUE.
      DO 3 K = 1,9
      H=A(K-1)-A(K)
    3 LMT=LMT .AND. (H .LE. 0 .OR. ABS(H) .LE. EPS*ABS(A(K)))
      IF(.NOT.LMT) GO TO 4

      DO 2 M = 1,9
      DO 2 K = 0,9-M
      IF(LEV(M)) THEN
       T(K,M)=W(M-1,1)*T(K+1,M-1)-W(M,1)*T(K,M-1)
      ELSEIF(LEV(K)) THEN
       T(K,M)=W(M-1,2)*T(K+1,M-1)-W(M,2)*T(K,M-1)
      ELSE
       T(K,M)=W(M-1,3)*T(K+1,M-1)-W(M,3)*T(K,M-1)
      ENDIF
    2 CONTINUE
      DFDX=T(0,9)
      RERR=0
      IF(DFDX .NE. 0) RERR=(DFDX-T(0,8))/DFDX
      DELTA=DEL
      RETURN
  101 FORMAT('FAILURE FOR X = ',1P,D15.8)
      END
