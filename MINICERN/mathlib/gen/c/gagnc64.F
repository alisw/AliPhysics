*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/04/01 15:02:05  mclareni
* Mathlib gen
*
*
#include "gen/pilot.h"
#if !defined(CERNLIB_DOUBLE)
      FUNCTION RGAGNC(A,X)
#endif
#if defined(CERNLIB_DOUBLE)
      FUNCTION DGAGNC(A,X)
#include "gen/imp64.inc"
#endif

C     Calculates the complementary incomplete gamma function G(A,X)
C     as defined in Ref. 1. Based on
C     1. W. Gautschi, ALGORITHM 542 Incomplete Gamma Functions,
C        ACM Trans. Math. Software 5 (1979) 482-489
C     2. W. Gautschi, A computational procedure for incomplete gamma
C        functions, ACM Trans. Math. Software 5 (1979) 466-481

      CHARACTER NAME*(*)
      CHARACTER*80 ERRTXT
#if !defined(CERNLIB_DOUBLE)
      PARAMETER (NAME = 'RGAGNC')
#endif
#if defined(CERNLIB_DOUBLE)
      PARAMETER (NAME = 'RGAGNC/DGAGNC')
#endif

      PARAMETER (EPS = 5D-14)
      PARAMETER (ALH = -0.69314 71805 59945 31D0)
      PARAMETER (Z1 = 1, HALF = Z1/2, QUAR = Z1/4)
      PARAMETER (C1 = 3*Z1/2, KMAX = 600, EPS1 = EPS/100)

      DIMENSION C(25)

      DATA C( 1) / 0.57721 56649 01532 86D0/
      DATA C( 2) /-0.65587 80715 20253 88D0/
      DATA C( 3) /-0.04200 26350 34095 24D0/
      DATA C( 4) / 0.16653 86113 82291 49D0/
      DATA C( 5) /-0.04219 77345 55544 34D0/
      DATA C( 6) /-0.00962 19715 27876 97D0/
      DATA C( 7) / 0.00721 89432 46663 10D0/
      DATA C( 8) /-0.00116 51675 91859 07D0/
      DATA C( 9) /-0.00021 52416 74114 95D0/
      DATA C(10) / 0.00012 80502 82388 12D0/
      DATA C(11) /-0.00002 01348 54780 79D0/
      DATA C(12) /-0.00000 12504 93482 14D0/
      DATA C(13) / 0.00000 11330 27231 98D0/
      DATA C(14) /-0.00000 02056 33841 70D0/
      DATA C(15) / 0.00000 00061 16095 10D0/
      DATA C(16) / 0.00000 00050 02007 64D0/
      DATA C(17) /-0.00000 00011 81274 57D0/
      DATA C(18) / 0.00000 00001 04342 67D0/
      DATA C(19) / 0.00000 00000 07782 26D0/
      DATA C(20) /-0.00000 00000 03696 81D0/
      DATA C(21) / 0.00000 00000 00510 04D0/
      DATA C(22) /-0.00000 00000 00020 58D0/
      DATA C(23) /-0.00000 00000 00005 35D0/
      DATA C(24) / 0.00000 00000 00001 23D0/
      DATA C(25) /-0.00000 00000 00000 12D0/

#if defined(CERNLIB_DOUBLE)
      GLGAMA(V)=DLGAMA(V)
#endif
#if !defined(CERNLIB_DOUBLE)
      GLGAMA(V)=ALGAMA(V)
#endif

      H=0
      IF(X .LT. 0) THEN
       WRITE(ERRTXT,101) X
       CALL MTLPRT(NAME,'C334.1',ERRTXT)
       GO TO 99
      ELSEIF(X .EQ. 0) THEN
       IF(A .LT. 0) THEN
        H=-1/A
       ELSEIF(A .EQ. 0) THEN
       CALL MTLPRT(NAME,'C334.2','ILLEGAL ARGUMENTS A = X = 0')
       ELSE
        H=1
       ENDIF
       GO TO 99
      ELSE
       ALX=LOG(X)
      ENDIF
      IF(X .LT. QUAR) THEN
       ALFA=ALH/ALX
      ELSE
       ALFA=X+QUAR
      ENDIF
      MA=HALF-A
      AEPS=A+MA

      IF(MA .GT. 0) THEN
       IF(AEPS .NE. 0) THEN
        ALGP1=GLGAMA(1+AEPS)-LOG(ABS(AEPS))
        IF(MA .NE. 1) ALGP1=ALGP1+GLGAMA(1-AEPS)-GLGAMA(MA-AEPS)
       ELSE
        ALGP1=0
       ENDIF
      ELSE
       ALGP1=GLGAMA(1+A)
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
    2  H=1-EXP(A*ALX-X+LOG(SUM)-ALGP1)
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
    4  IF(A .LE. 0) THEN
        H=SUM/(X+1-A)
       ELSE
        H=EXP(A*ALX-X+LOG(A*SUM/(X+1-A))-ALGP1)
       ENDIF
      ELSE
       AE=A
       IF(A .LT. HALF) THEN
        IF(A .LT. -HALF) AE=AEPS
        SUM=C(25)
        DO 12 K = 24,1,-1
   12   SUM=AE*SUM+C(K)
        GA=-SUM/(1+AE*SUM)
        Y=AE*ALX
        IF(ABS(Y) .GE. 1) THEN
         U=GA-(EXP(Y)-1)/AE
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
        U=EXP(GLGAMA(A))-X**A/A
       ENDIF
       P=AE*X
       Q=AE+1
       R=AE+3
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
   10  H=U+SUM*X**(AE+1)/(AE+1)
       IF(A .LT. -HALF) THEN
        H=H*EXP(X-AE*ALX)
        DO 13 J = 1,MA
   13   H=(1-X*H)/(J-AE)
       ELSEIF(A .LE. 0) THEN
        H=H*EXP(X-A*ALX)
       ELSE
        H=A*H*EXP(-ALGP1)
       ENDIF
      ENDIF
#if defined(CERNLIB_DOUBLE)
   99 DGAGNC=H
#endif
#if !defined(CERNLIB_DOUBLE)
   99 RGAGNC=H
#endif
      RETURN

   98 WRITE(ERRTXT,103) A,X
      CALL MTLPRT(NAME,'C334.3',ERRTXT)
      GO TO 99
  101 FORMAT('ILLEGAL ARGUMENT  X = ',1P,E15.6,' < 0')
  103 FORMAT('PROBLEMS WITH CONVERGENCE, A = ',1P,E15.8,'  X = ',E15.6)
      END
