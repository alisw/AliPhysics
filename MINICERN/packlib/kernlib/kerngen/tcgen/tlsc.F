*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/02/15 17:49:53  mclareni
* Kernlib
*
*
#include "kerngen/pilot.h"
      SUBROUTINE TLSC (A,B,AUX,IPIV,EPS,X)
C
C CERN PROGLIB# E230    TLSC            .VERSION KERNFOR  2.06  740511
C ORIG. 11/05/74 WH+WM
C
C.  SUBROUTINE TLSC           CONSTRAINED L.S.              HART/MATT
C.
C.        M1 PARAMETERS ARE ELIMINATED USING THE CONSTRAINTS AND
C.        THE REMAINING DERIVATIVE PART OF MTX A TRIANGULARISED BY
C.        HOUSEHOLDER ROTATIONS. X IS FOUND BY BACK SUBSTITUTION.
C.  ARGUMENTS
C.        A    M BY N CONSTRAINT/COEFFICIENT MTX (DESTROYED)
C.                    (CONSTRAINTS FIRST).
C.        B    M BY L R.H.S. MTX (DESTROYED)
C.        AUX  AUXILIARY STORAGE ARRAY OF DIM=MAX(N,L)+N
C.             ON RETURN THE FIRST L LOCS CONTAIN THE RESULTING
C.             SUM OF SQUARES
C.        IPIV INTEGER VECTOR (DIM=N) WHICH RETURNS THE COLUMN
C.             INTERCHANGES IN MTX A
C.        EPS  INPUT PARAM SPECIFYING A TOLERANCE FOR PIVOTING.
C.             NO EXCHANGE OF COL(I) UNLESS EPS*PIV(I) .GT. PIV(1).
C.        X    N BY L SOLUTION MATRIX.
C.        THE FOLLOWING ARE TRANSMITTED THRO COMMON /TLSDIM/  ...
C.        M1          NUMBER OF CONSTRAINT EQUATIONS.
C.        M           TOTAL NUMBER OF EQUATIONS.
C.        N           NUMBER OF UNKNOWNS.
C.        L           NUMBER OF SYSTEMS TO BE SOLVED.
C.        IER  OUTPUT PARAM GIVING RANK OF MTX A
C.             IF VALUE IS -1001 NO SOLUTION CAN BE FOUND.
C.  REMARKS
C.        ONLY PIVOT IF NECESSARY
C.
C.-------------------------------------------------------------------
C
      COMMON /TLSDIM/ M1,M,N,L,IER
      COMMON /SLATE/ BETA,H,I,IB,IB1,ID,ID1,IEND,II,IST,J,JA,JB,JK
     +              ,JST,K,KPIV,KR,KST,KT,K1,LV,MR,M11,NK,NR,PIV,PIVT
     +              ,SIG,DUM(11)
      DIMENSION      A(*), AUX(*), B(*), IPIV(*), X(*)
C
C--     ERROR TEST
      IF (N.GT.M.OR.M1.GT.N)  GO TO 90
C--      CALCULATION OF INITIAL VECTORS.
C--
      K1  = MAX (N,L)
      IER = 1
      DO           5         K=1,N
    5 IPIV(K) = K
C
C--      PERFORM A DECOMPOSITION LOOP.
C
      IST = - N
      JB  = 1 - L
      M11 = M1 + 1
      MR  = M1
      DO           50        K=1,N
      IF     (K.GT.M1)                 MR = M
      IST = IST + N + 1
      JB  = JB + L
      LV  = MR - K + 1
C
C--      PIVOT ONLY IF ABSOLUTLY NECESSARY.
C
      PIV = 0.
      ID  = IST - N
      DO           20        J=K,N
      IF     (K.EQ.1 .OR. K.EQ.M11)    GO TO     10
      PIVT = AUX(J) - A(ID)*A(ID)
      GO TO      15
C
   10 I = ID + N
      IF (LV .EQ. 1) GO TO 12
      CALL TLSMSQ (A(I),N,LV,PIVT)
      GO TO 15
   12 PIVT= A(I)*A(I)
C
   15 AUX(J) = PIVT
      ID = ID + 1
      IF     (PIVT*EPS.LE.PIV)         GO TO     20
      PIV  = PIVT
      KPIV = J
   20 CONTINUE
C
      I = KPIV - K
      IF     (I.LE.0)                  GO TO     25
      H = AUX(K)
      AUX(K) = AUX(KPIV)
      AUX(KPIV) = H
      ID = IST + I
      NR = M - K + 1
      CALL TLSWOP (A(IST),A(ID),N,NR)
C
C--      COMPUTATION OF TRANSFORMATION PARAMETER SIG.
C--      GENERATION OF VECTOR UK IN K-TH COLUMN OF MATRIX A AND OF
C--      PARAMETER BETA.
C
   25 CALL TLUK (A(IST),N,LV,SIG,BETA)
      IF     (LV.EQ.0)                 GO TO     90
C
      J = K1 + K
      AUX(J)=-SIG
      IF     (K.GE.N)                  GO TO     30
C
C--      TRANSFORMATION OF MATRIX A.
C
      NK = N - K
      IF (LV.EQ.1)                     GO TO     27
      CALL TLSTEP (A(IST),A(IST+1),N,N,LV,NK,BETA)
      GO TO        30
   27 DO           28        J=1,NK
      JST = IST + J
   28 A(JST) = A(JST)*(1.-BETA*A(IST)**2)
C
C--      TRANSFORMATION OF RIGHT HAND SIDE MATRIX B.
C
   30 IB = (K-1) * L + 1
      IF (LV.EQ.1)                     GO TO     32
      CALL TLSTEP (A(IST),B(IB),N,L,LV,L,BETA)
      GO TO      34
   32 DO           33        J=1,L
      JST = IB + J - 1
   33 B(JST) = B(JST)*(1.-BETA*A(IST)**2)
   34 IPIV(KPIV) = IPIV(K)
      IPIV(K) = KPIV
      IF     (K.GT.M1)                 GO TO     50
C
C--      REDUCE MATRIX A FOR ONE CONSTRAINT PARAMETER.
C
      DO           45        I=M11,M
      ID1 = IST + (I-K)*N
      IF     (A(ID1).EQ.0)             GO TO     45
      H  = - A(ID1)/SIG
      A(ID1) = H
      ID1 = ID1 + 1
      ID = IST + 1
C
      DO           35        J=1,NK
      A(ID1) = A(ID1) - H*A(ID)
      ID1 = ID1 + 1
   35 ID  = ID + 1
C
      IB1 = 1 + (I-1)*L
      IB = JB
      DO           40        J=1,L
      B(IB1) = B(IB1) - H*B(IB)
      IB1 = IB1 + 1
   40 IB  = IB + 1
C
   45 CONTINUE
   50 CONTINUE
C
C--      BACK SUBSTITUTION AND BACK INTERCHANGE.
C
      IER = N * IER
      KT  = N
      JK  = (N-1) * L
      K   = K1 + N
      PIV = 1./AUX(K)
      DO           55        K=1,L
      JK = JK + 1
   55 X(JK) = PIV * B(JK)
C
      KR = N - 1
      IF     (KR.LE.0)                 GO TO     70
C
      JST = KR * (N+1) + 2
      DO           65        J=1,KR
      JST = JST - N - 1
      IEND= (KR-J+1) * N
      K   = K1 + KR - J + 1
      PIV = 1./AUX(K)
      KST = K-K1
      ID  = IPIV(KST)-KST
      KST = (KR-J) * L
C
      DO           65        K=1,L
      KST = KST + 1
      H=B(KST)
C
      II = KST
      DO           60        I=JST,IEND
      II = II + L
   60 H  = H - A(I) * X(II)
C
      II = KST + ID *L
      X(KST) = X(II)
      X(II)  = PIV * H
   65 CONTINUE
C
C--      COMPUTATION OF LEAST SQUARES.
C
   70 IST = KT*L
      DO           80        J=1,L
      IST = IST + 1
      H   = 0.
      JA  = IST
      IF     (M.LE.KT)                 GO TO     80
      NR  = M - KT
      IF     (NR.EQ.1)                 GO TO     75
      CALL TLSMSQ (B(IST),L,NR,H)
      GO TO      80
C
   75 H = B(IST)*B(IST)
C
   80 AUX(J) = H
      RETURN
C--
C--      ERROR RETURN IN CASE OF A ZERO COLUMN IN MATRIX A.
C
   90 IER = -1001
      RETURN
      END
