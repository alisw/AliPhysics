*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/04/01 15:01:52  mclareni
* Mathlib gen
*
*
#include "gen/pilot.h"
#if defined(CERNLIB_DOUBLE)
      SUBROUTINE DSNLEQ(N,X,F,FTOL,XTOL,MAXF,IPRT,INFO,SUB,W)
#include "gen/imp64.inc"
#endif
#if !defined(CERNLIB_DOUBLE)
      SUBROUTINE  RSNLEQ(N,X,F,FTOL,XTOL,MAXF,IPRT,INFO,SUB,W)
#endif

C     Based on   J.J. More  and  M.Y. Cosnard
C
C       ALGORITHM 554 BRENTM, A Fortran Subroutine for the
C       Numerical Solution of Systems of Nonlinear Equations [C5]
C
C     ACM Trans. Math. Software 6 (1980) 240-251.

      DIMENSION X(N),F(N),W(N,*),MPT(288)
      LOGICAL LCV

      PARAMETER (Z1 = 1, SCALE = 10, P05 = 5*Z1/100)
C**** EPS = SQRT(SMALLEST FP.NUMBER)
C     EPS = 1 / SQRT( 16D0**13 )
#if !defined(CERNLIB_DOUBLE)
      PARAMETER (EPS =  0.84293 69702 17878 97282 52636 392E-07)
#endif
#if defined(CERNLIB_IBM)
      PARAMETER (EPS =  0.14901 16119 38476 562D-07)
#endif
#if defined(CERNLIB_VAX)
      PARAMETER (EPS =  0.37252 90298 46191 40625D-08)
#endif
#if (defined(CERNLIB_UNIX))&&(defined(CERNLIB_DOUBLE))
      PARAMETER (EPS =  0.14901 16119 38476 600D-07)
#endif
      DATA (MPT(I),I=1,288)
     1/1* 1,1* 2,3* 3,3* 4,4* 5,4* 6,4* 7,4* 8,5* 9,5*10,5*11,5*12,
     2 5*13,5*14,6*15,6*16,5*17,6*18,6*19,6*20,7*21,6*22,6*23,7*24,
     3 6*25,7*26,6*27,7*28,7*29,7*30,7*31,7*32,7*33,7*34,7*35,7*36,
     4 8*37,7*38,7*39,8*40,7*41,8*42,7*43,8*44,8*45,7*46,8*47,8*48/

      INFO=0
      IF(N .LE. 0 .OR. FTOL .LE. 0 .OR. XTOL .LE. 0) RETURN
C
C     Find optimal MOPT for iterative refinement
C
      IF(N .LE. 288) THEN
       MOPT=MPT(N)
      ELSE
       H=0
       DO 1 I = 49,N
       TEMP=LOG(I+Z1)/(N+2*I+1)
       IF(TEMP .LT. H) THEN
        MOPT=I-1
        GO TO 2
       ENDIF
    1  H=TEMP
      ENDIF

    2 IFLAG=0
      NUMF=0
      NFCALL=0

      NIER6=-1
      NIER7=-1
      NIER8=0
      FNORM=0
      DIFIT=0
      XNORM=0
      DO 10 I = 1,N
   10 XNORM=MAX(XNORM,ABS(X(I)))
      DELTA=SCALE*XNORM
      IF(XNORM .EQ. 0) DELTA=SCALE

   20 IF(IPRT .NE. 0) WRITE(6,'(1X,I5,D25.14)') (I,X(I),I=1,N)

      NSING=N
      FNORM1=FNORM
      DIFIT1=DIFIT
      FNORM=0
C
C     Compute step H for the divided difference which approximates
C     the K-th row of the Jacobian matrix
C
      H=EPS*XNORM
      IF(H .EQ. 0) H=EPS
      DO 40 J = 1,N
      DO 30 I = 1,N
   30 W(I,J+3)=0
      W(J,J+3)=H
   40 W(J,2)=X(J)
C
C     Enter a subiteration
C
      DO 150 K = 1,N
      IFLAG=K
      CALL SUB(N,W(1,2),F,IFLAG)
      FKY=F(K)
      NFCALL=NFCALL+1
      NUMF=NFCALL/N
      IF(IFLAG .LT. 0) GO TO 230
      FNORM=MAX(FNORM,ABS(FKY))
C
C     Compute the K-th row of the Jacobian matrix
C
      DO 60 J = K,N
      DO 50 I = 1,N
   50 W(I,3)=W(I,2)+W(I,J+3)
      CALL SUB(N,W(1,3),F,IFLAG)
      FKZ=F(K)
      NFCALL=NFCALL+1
      NUMF=NFCALL/N
      IF(IFLAG .LT. 0) GO TO 230
   60 W(J,1)=FKZ-FKY
      F(K)=FKY
C
C     Compute the Householder transformation to reduce the K-th row
C     of the Jacobian matrix to a multiple of the K-th unit vector
C
      ETA=0
      DO 70 I = K,N
   70 ETA=MAX(ETA,ABS(W(I,1)))
      IF(ETA .EQ. 0) GO TO 150
      NSING=NSING-1
      SKNORM=0
      DO 80 I = K,N
      W(I,1)=W(I,1)/ETA
   80 SKNORM=SKNORM+W(I,1)**2
      SKNORM=SQRT(SKNORM)
      IF(W(K,1) .LT. 0) SKNORM=-SKNORM
      W(K,1)=W(K,1)+SKNORM
C
C     Apply the transformation
C
      DO 90 I = 1,N
   90 W(I,3)=0
      DO 100 J = K,N
      DO 100 I = 1,N
  100 W(I,3)=W(I,3)+W(J,1)*W(I,J+3)
      DO 120 J = K,N
      TEMP=W(J,1)/(SKNORM*W(K,1))
      DO 120 I = 1,N
  120 W(I,J+3)=W(I,J+3)-TEMP*W(I,3)
C
C     Compute the subiterate
C
      W(K,1)=SKNORM*ETA
      TEMP=FKY/W(K,1)
      IF(H*ABS(TEMP) .GT. DELTA) TEMP=SIGN(DELTA/H,TEMP)
      DO 140 I = 1,N
  140 W(I,2)=W(I,2)+TEMP*W(I,K+3)
  150 CONTINUE
C
C     Compute the norms of the iterate and correction vector
C
      XNORM=0
      DIFIT=0
      DO 160 I = 1,N
      XNORM=MAX(XNORM,ABS(W(I,2)))
      DIFIT=MAX(DIFIT,ABS(X(I)-W(I,2)))
  160 X(I)=W(I,2)
C
C     Update the bound on the correction vector
C
      DELTA=MAX(DELTA,SCALE*XNORM)
C
C     Determine the progress of the iteration
C
      LCV=FNORM .LT. FNORM1 .AND. DIFIT .LT. DIFIT1 .AND. NSING .EQ. 0
      NIER6=NIER6+1
      NIER7=NIER7+1
      NIER8=NIER8+1
      IF(LCV) NIER6=0
      IF(FNORM .LT. FNORM1 .OR. DIFIT .LT. DIFIT1) NIER7=0
      IF(DIFIT .GT. EPS*XNORM) NIER8=0
C
C     Tests for convergence
C
      IF(FNORM .LE. FTOL) INFO=1
      IF(DIFIT .LE. XTOL*XNORM .AND. LCV) INFO=2
      IF(FNORM .LE. FTOL .AND. INFO .EQ. 2) INFO=3
      IF(INFO .NE. 0) GO TO 230
C
C     Tests for termination
C
      IF(NUMF .GE. MAXF) INFO=4
      IF(NSING .EQ. N) INFO=5
      IF(NIER6 .EQ. 5) INFO=6
      IF(NIER7 .EQ. 3) INFO=7
      IF(NIER8 .EQ. 4) INFO=8
      IF(INFO .NE. 0) GO TO 230
      IF(.NOT.LCV .OR. DIFIT .GT. P05*XNORM) GO TO 20
C
C     Iterative refinement  (if the iteration is converging)
C
      DO 210 M = 2,MOPT
      FNORM1=FNORM
      FNORM=0
      DO 190 K = 1,N
      IFLAG=K
      CALL SUB(N,W(1,2),F,IFLAG)
      FKY=F(K)
      NFCALL=NFCALL+1
      NUMF=NFCALL/N
      IF(IFLAG .LT. 0) GO TO 230
      FNORM=MAX(FNORM,ABS(FKY))
C
C     Iterative refinement is terminated if it does not give a
C     reduction on residuals
C
      IF(FNORM .GE. FNORM1) THEN
       FNORM=FNORM1
       GO TO 20
      ENDIF
      TEMP=FKY/W(K,1)
      DO 180 I = 1,N
  180 W(I,2)=W(I,2)+TEMP*W(I,K+3)
  190 CONTINUE
C
C     Compute the norms of the iterate and correction vector
C
      XNORM=0
      DIFIT=0
      DO 200 I = 1,N
      XNORM=MAX(XNORM,ABS(W(I,2)))
      DIFIT=MAX(DIFIT,ABS(X(I)-W(I,2)))
  200 X(I)=W(I,2)
C
C     Stopping criteria for iterative refinement
C
      IF(FNORM .LE. FTOL) INFO=1
      IF(DIFIT .LE. XTOL*XNORM) INFO=2
      IF(FNORM .LE. FTOL .AND. INFO .EQ. 2) INFO=3
      IF(NUMF .GE. MAXF .AND. INFO .EQ. 0) INFO=4
      IF(INFO .NE. 0) GO TO 230
  210 CONTINUE
      GO TO 20

  230 IF(IFLAG .LT. 0) INFO=IFLAG
      RETURN
      END
