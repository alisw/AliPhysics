*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/02/15 17:49:53  mclareni
* Kernlib
*
*
#include "kerngen/pilot.h"
      SUBROUTINE TLS (A,B,AUX,IPIV,EPS,X)
C
C CERN PROGLIB# E230    TLS             .VERSION KERNFOR  2.06  740511
C ORIG. 11/05/74 W.HART+W.MATT
C
C.
C.  SUBROUTINE TLS            LINEAR LEAST SQUARES          HART/MATT
C.        HOUSEHOLDER ROTATIONS ARE USED TO ROTATE A TO TRIANGULAR
C.        FORM THEN X OBTAINED THRO BACK SUBSTITUTION.
C.        IDENTICAL IN OPERATION TO LIB ROUTINE HLS BUT WITH TC ARRAY
C.        DEFINITIONS. SEE HLS WRITE UP (D520) FOR FURTHER DETAILS.
C.  ARGUMENTS
C.        A    M BY N CONSTRAINT/COEFFICIENT MTX (DESTROYED),
C.             CONSTRAINTS FIRST.
C.        B    M BY L R.H.S. MTX (DESTROYED)
C.        AUX  AUXILIARY STORAGE ARRAY OF DIM=MAX(N,L)+N
C.             ON RETURN THE FIRST L LOCS CONTAIN THE RESULTING
C.             SUM OF SQUARES
C.        IPIV INTEGER VECTOR (DIM=N) WHICH RETURNS THE COLUMN
C.             INTERCHANGES IN MTX A
C.        EPS  INPUT PARAM SPECIFYING A TOLERANCE FOR DETERMINING
C.             THE RANK OF MTX A
C.        X    N BY L SOLUTION MTX
C.        THE FOLLOWING ARE TRANSMITTED THRO COMMON /TLSDIM/  ...
C.        M1   NUMBER OF CONSTRAINT EQUATIONS = 0
C.        M    TOTAL NUMBER OF EQUATIONS
C.        N    NUMBER OF UNKNOWNS
C.        L    NUMBER OF SYSTEMS TO BE SOLVED
C.        IER  OUTPUT PARAM GIVING RANK OF MTX A
C.
C.-------------------------------------------------------------------
C
      COMMON /TLSDIM/ M1,M,N,L,IER
      COMMON /SLATE/ BETA,F,G,H,I,IB,ID,IEND,II,IST,J,JA,JB,JD,JK,JL
     +              ,TOL,JX,J1,K,KPIV,KR,KST,KT,K1,K2,K3,LV,NK,NR,N1
     +              ,PIV,PIVT,SIG,JST,DUM(5)
      DIMENSION A(*),B(*),X(*),IPIV(*),AUX(*)
C
C     ERROR TEST
      IF(M-N)31,1,1
C
C     CALCN OF INITIAL VECTORS S(K),T(K) IN LOCS AUX(1),AUX(K1)
    1 PIV=0.
      K1 = MAX (N,L)
      IST = 0
      K2 = K1 + 1
      DO 4 K=1,N
      IPIV(K)=K
      IST = IST + 1
      IF (M .EQ. 1) GO TO 40
      CALL TLPIV(A(IST),B,N,L,M,H,G)
      GO TO 41
   40 G=A(IST)*B(1)
      H=A(IST)*A(IST)
   41 AUX(K)=H
      AUX(K2) = G
      PIVT = G*G/H
      IF(PIVT-PIV)4,4,3
    3 PIV = PIVT
      KPIV=K
    4 K2 = K2 + 1
      IF (M .EQ. 1) GO TO 2
      CALL TLSMSQ(B(1),L,M,F)
      GO TO 7
    2 F=B(1)*B(1)
C
C     ERROR TEST
    7 IF(PIV)31,31,5
C
C     DEFINE TOLERANCE FOR CHECKING RANK OF A
    5 TOL = EPS*EPS
      IER = 1
C
C
C     DECOMPOSITION LOOP
    9 IST = -N
      JB = 1 - L
      DO 21 K=1,N
      IF(EPS.EQ.0.) GO TO 12
      IF(EPS.GT.0.) GO TO 11
      IF(PIV.GT.TOL) GO TO 12
      GO TO 34
   11 IF(F.GT.TOL*FLOAT(M-K+1))  GO TO 12
   34 IF(K.NE.1) GO TO 32
      TOL = 0.
      IER = -1
   12 F = F - PIV
      IST = IST + N + 1
      JB = JB + L
      LV = M-K+1
      I=KPIV-K
      IF(I)8,8,6
C
C     INTERCHANGE K-TH COLUMN OF A WITH KPIV-TH IN CASE KPIV.GT.K
    6 H=AUX(K)
      AUX(K)=AUX(KPIV)
      AUX(KPIV)=H
      K2 = K1 + K
      K3 = K1 + KPIV
      G = AUX(K2)
      AUX(K2) = AUX(K3)
      AUX(K3) = G
      JA = IST
      JD = IST + I
      NR = M - K + 1
      CALL TLSWOP(A(IST),A(JD),N,NR)
C
C     COMPUTATION OF PARAMETER SIG
C     GENERATION OF VECTOR UK IN K-TH COLUMN OF MATRIX A AND OF
C     PARAMETER BETA
    8 CALL TLUK(A(IST),N,LV,SIG,BETA)
C
      J = K1 + K
      AUX(J)=-SIG
C
C     SAVE INTERCHANGE INFORMATION
   13 IPIV(KPIV)=IPIV(K)
      IPIV(K)=KPIV
      IF(K-N)14,19,19
C
C     TRANSFORMATION OF MATRIX A
   14 NK = N - K
      CALL TLSTEP(A(IST),A(IST+1),N,N,LV,NK,BETA)
C
C     TRANSFORMATION OF RIGHT HAND SIDE MATRIX B
   19 IB = (K-1) * L +1
      CALL TLSTEP(A(IST),B(IB),N,L,LV,L,BETA)
      IF(K-N)10,21,21
C
C     UPDATING OF S(K),T(K) ELEMENTS STORED IN AUX
   10 PIV = 0.
      KPIV = K + 1
      J1 = KPIV
      ID = IST
      K2 = K1 + J1
      DO 18 J=J1,N
      ID = ID+1
      H=AUX(J) - A(ID) * A(ID)
      AUX(J)=H
      G = AUX(K2) - A(ID) * B(JB)
      AUX(K2) = G
      PIVT = G*G/H
      IF(PIVT-PIV)18,18,17
   17 PIV = PIVT
      KPIV=J
   18 K2 = K2 + 1
   21 CONTINUE
      GO TO 20
C     END OF DECOMPOSITION LOOP
C
C
C     RANK OF MATRIX LESS THAN N, ZERO X,S AND BACK SUBSTITUTE
   32 KR = K - 1
      KT = KR
      IER = KR
      JK = KR * L + 1
      JL = N * L
      DO 15 JX=JK,JL
   15 X(JX) = 0.
      GO TO 16
C
C     BACK SUBSTITUTION AND BACK INTERCHANGE
   20 IER = N * IER
      KT = N
      JK = (N-1) * L
      K = K1 + N
      PIV=1./AUX(K)
      DO 22 K=1,L
      JK = JK + 1
   22 X(JK) = PIV * B(JK)
      KR = N - 1
C
   16 IF(KR)26,26,23
   23 JST = KR * (N+1) + 2
      DO 25 J=1,KR
      JST = JST - N - 1
      IEND = (KR-J+1) * N
      K = K1 + KR - J + 1
      PIV=1./AUX(K)
      KST=K-K1
      ID=IPIV(KST)-KST
      KST = (KR-J) * L
      DO 25 K=1,L
      KST = KST + 1
      H=B(KST)
      II = KST
      DO 24 I = JST,IEND
      II = II + L
   24 H = H - A(I) * X(II)
      II = KST + ID * L
      X(KST) = X(II)
      X(II)=PIV*H
   25 CONTINUE
C
C     COMPUTATION OF LEAST SQUARES
   26 IST = KT * L
      N1 = KT + 1
      DO 29 J=1,L
      IST = IST + 1
      H=0.
      JA = IST
      IF(M-KT)29,29,27
   27 NR = M - N1 + 1
      IF(NR.EQ.1) GO TO 28
      CALL TLSMSQ(B(IST),L,NR,H)
      GO TO 29
   28 H = B(IST)*B(IST)
C
   29 AUX(J)=H
      RETURN
C
C     ERROR RETURN IN CASE OF ZERO-MATRIX A OR M.LT.N
   31 IER = -1001
      RETURN
      END
