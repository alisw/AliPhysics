*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/04/01 15:02:27  mclareni
* Mathlib gen
*
*
#include "gen/pilot.h"
      SUBROUTINE CHEB(M,N,MDIM,NDIM,A,B,TOL,RELERR,X,RANK,RESMAX,
     +   ITER,OCODE)
C         LINEAR CHEBYSHEV FIT ROUTINE FROM ACM TOMS 1, 266(1975)
C         BY BARRODALE AND PHILLIPS
      DIMENSION A(NDIM,MDIM), B(MDIM), X(NDIM)
      INTEGER PROW,PCOL,RANK,RANKP1,OCODE
#include "e221prec.inc"
      MP1 = M + 1
      NP1 = N + 1
      NP2 = N + 2
      NP3 = N + 3
      NP1MR = 1
      RANK = N
      RELTMP = RELERR
      RELERR = 0.
      DO 10 J= 1,M
        A(NP1,J) = 1.
        A(NP2,J) = -B(J)
        A(NP3,J) = N+J
   10 CONTINUE
      A(NP1,MP1) = 0.
      ITER = 0
      OCODE = 1
      DO 20 I= 1, N
        X(I)=0.
        A(I,MP1)=I
   20 CONTINUE
C LEVEL 1.
      LEV = 1
      K = 0
   30 K = K + 1
      KP1 = K + 1
      NP1MK = NP1 - K
      MODE = 0
      CALL VFILL(B(K),M-K+1,1.)
C DETERMINE THE VECTOR TO ENTER THE BASIS.
   50 D = -BIG
      DO 60 J= K, M
        IF (B(J) .EQ. 0.)  GOTO 60
        DD = ABS(A(NP2,J))
        IF (DD .LE. D)  GOTO 60
        PCOL = J
        D = DD
   60 CONTINUE
      IF (K .GT. 1) GOTO 70
C TEST FOR ZERO RIGHT-HAND SIDE
      IF (D .GT. TOL)  GOTO 70
      RESMAX = 0.
      MODE = 2
      GOTO 380
C DETERMINE THE VECTOR TO LEAVE THE BASIS.
   70 D = TOL
      DO 80 I= 1, NP1MK
        DD = ABS(A(I,PCOL))
        IF (DD .LE. D) GOTO 80
        PROW = I
        D = DD
   80 CONTINUE
      IF (D .GT. TOL) GOTO 330
C CHECK FOR LINEAR DEPENDENCE IN LEVEL 1.
      B(PCOL) = 0.
      IF (MODE .EQ. 1)  GOTO 50
      DO 100 J= K, M
        IF (B(J) .EQ. 0.)  GOTO 100
        DO 90 I= 1, NP1MK
          IF (ABS(A(I,J)) .LE. TOL)  GOTO 90
          MODE = 1
          GOTO 50
   90   CONTINUE
  100 CONTINUE
      RANK = K - 1
      NP1MR = NP1 - RANK
      OCODE = 0
      GOTO 160
  110 IF (PCOL .EQ. K)  GOTO 130
C INTERCHANGE COLUMNS IN LEVEL 1.
      DO 120 I= 1, NP3
        D = A(I,PCOL)
        A(I,PCOL) = A(I,K)
        A(I,K) = D
  120 CONTINUE
  130 IF (PROW .EQ. NP1MK)  GOTO 150
C INTERCHANGE ROWS IN LEVEL 1
      DO 140 J= 1, MP1
        D = A(PROW,J)
        A(PROW,J) = A(NP1MK,J)
        A(NP1MK,J) = D
  140 CONTINUE
  150 IF (K.LT.N)  GOTO 30
  160 IF (RANK.EQ.M) GOTO 380
      RANKP1 = RANK + 1
C LEVEL 2.
      LEV = 2
C DETERMINE THE VECTOR TO ENTER THE BASIS.
      D = TOL
      DO 170 J= RANKP1, M
        DD = ABS(A(NP2,J))
        IF (DD .LE. D)  GOTO 170
        PCOL = J
        D = DD
  170 CONTINUE
C COMPARE CHEBYSHEV ERROR WITH TOL.
      IF (D.GT.TOL)  GOTO 180
      RESMAX = 0.
      MODE = 3
      GOTO 380
C
  180 IF (A(NP2,PCOL) .LT. -TOL)  GOTO 200
      A(NP1,PCOL) = 2. - A(NP1,PCOL)
      DO 190 I= NP1MR, NP3
        IF (I.EQ.NP1) GOTO 190
        A(I,PCOL) = -A(I,PCOL)
  190 CONTINUE
C ARRANGE FOR ALL ENTRIES IN PIVOT COL (EXC. PIVOT) TO BE NEGATIVE
  200 DO 220 I= NP1MR, N
        IF (A(I,PCOL) .LT. TOL) GOTO 220
        DO 210 J= 1, M
          A(NP1,J) = A(NP1,J) + 2.*A(I,J)
          A(I,J) = -A(I,J)
  210   CONTINUE
        A(I,MP1) = -A(I,MP1)
  220 CONTINUE
      PROW = NP1
      GOTO 330
C
  230 IF (RANKP1 .EQ. M)  GOTO 380
      IF (PCOL .EQ. M)  GOTO 250
C INTERCHANGE COLUMNS IN LEVEL 2.
      DO 240 I= NP1MR,NP3
        D = A(I,PCOL)
        A(I,PCOL) =A(I,M)
        A(I,M) = D
  240 CONTINUE
  250 MM1 = M-1
C LEVEL 3.
      LEV = 3
C DETERMINE THE VECTOR TO ENTER THE BASIS.
  260 D= -TOL
      VAL = 2. * A(NP2,M)
      DO 280 J= RANKP1,MM1
        IF (A(NP2,J) .GE. D)  GOTO 270
        PCOL = J
        D = A(NP2,J)
        MODE = 0
        GOTO 280
  270   DD = VAL - A(NP2,J)
        IF (DD .GE. D)  GOTO 280
        MODE = 1
        PCOL = J
        D = DD
  280 CONTINUE
      IF (D .GE. -TOL)  GOTO 380
      DD = -D/A(NP2,M)
      IF (DD .GE. RELTMP)  GOTO 290
      RELERR = DD
      MODE = 4
      GOTO 380
C
  290 IF (MODE .EQ. 0)  GOTO 310
      DO 300 I= NP1MR,NP1
        A(I,PCOL) = 2.*A(I,M) - A(I,PCOL)
  300 CONTINUE
      A(NP2,PCOL) = D
      A(NP3,PCOL) = -A(NP3,PCOL)
C DETERMINE THE VECTOR TO LEAVE THE BASIS.
  310 D = BIG
      DO 320 I= NP1MR,NP1
        IF (A(I,PCOL) .LE. TOL)  GOTO 320
        DD = A(I,M)/A(I,PCOL)
        IF (DD .GE. D) GOTO 320
        PROW = I
        D = DD
  320 CONTINUE
      IF (D .LT. BIG)  GOTO 330
      OCODE = 2
      GOTO 380
C PIVOT ON A(PROW,PCOL)
  330 PIVOT = A(PROW,PCOL)
      DO 340 J= 1, M
        A(PROW,J) = A(PROW,J)/PIVOT
  340 CONTINUE
      DO 360 J= 1, M
        IF (J.EQ.PCOL) GOTO 360
        D = A(PROW,J)
        DO 350 I= NP1MR,NP2
          IF (I.EQ.PROW)  GOTO 350
          A(I,J) = A(I,J) - D*A(I,PCOL)
  350   CONTINUE
  360 CONTINUE
      TPIVOT = -PIVOT
      DO 370 I= NP1MR, NP2
        A(I,PCOL) = A(I,PCOL)/TPIVOT
  370 CONTINUE
      A(PROW,PCOL) = 1./PIVOT
      D = A(PROW,MP1)
      A(PROW,MP1) = A(NP3,PCOL)
      A(NP3,PCOL) = D
      ITER = ITER + 1
      GOTO (110, 230, 260), LEV
C
C PREPARE OUTPUT
  380 DO 390 J= 1, M
        B(J) = 0.
  390 CONTINUE
      IF (MODE .EQ. 2)  GOTO 450
      DO 400 J= 1, RANK
        K= A(NP3,J)
        X(K) = A(NP2,J)
  400 CONTINUE
      IF (MODE.EQ.3 .OR. RANK.EQ.M)  GOTO 450
      DO 410 I= NP1MR,NP1
        K = ABS(A(I,MP1)) - N
        B(K) = A(NP2,M) * SIGN(1., A(I,MP1))
  410 CONTINUE
      IF (RANKP1 .EQ. M)  GOTO 430
      DO 420 J= RANKP1, MM1
        K = ABS(A(NP3,J)) - N
        B(K) = (A(NP2,M)-A(NP2,J)) * SIGN(1., A(NP3,J))
  420 CONTINUE
C TEST FOR NON-UNIQUE SOLUTION
  430 DO 440 I= NP1MR, NP1
        IF (ABS(A(I,M)) .GT. TOL)  GOTO 440
        OCODE = 0
        GOTO 450
  440 CONTINUE
  450 IF (MODE.NE.2 .AND. MODE.NE.3) RESMAX = A(NP2,M)
      IF (RANK.EQ.M) RESMAX=0.
      IF (MODE.EQ.4) RESMAX = RESMAX-D
      RETURN
      END
