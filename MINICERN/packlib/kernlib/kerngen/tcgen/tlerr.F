*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/02/15 17:49:53  mclareni
* Kernlib
*
*
#include "kerngen/pilot.h"
      SUBROUTINE TLERR (A,X,AUX,IPIV)
C
C CERN PROGLIB# E230    TLERR           .VERSION KERNFOR  2.06  740511
C ORIG. 11/05/74 WH+WM
C
C.  SUBROUTINE TLERR          L.S. ERROR MATRIX             HART/MATT
C.
C.        CALCULATES INVERSE OF (A-TRANSPOSED*A) DIRECTLY FROM THE
C.        TRIANGULARISED TRANSFORM OF A.
C.  ARGUMENTS
C.        A,X,AUX,IPIV,(M1,M,N,L,IER) DEFINED AS FOR TLSC WITH X THE
C.        COVARIANCE MATRIX. MATRIX X MAY OVERWRITE MATRIX A.
C.  REMARKS
C.        CONSTRAINED AND UNCONSTRAINED VERSION COMBINED
C.
C.-------------------------------------------------------------------
C
      COMMON /TLSDIM/ M1,M,N,L,IER
      COMMON /SLATE/ BETA,H,I,IA,ID,IEND,II,IK,IL,IST,J,JA,JD,JK,JL
     +              ,K,KN,KS,K1,LV,NK,PIV,SIG,DUM(17)
      DIMENSION      A(*), AUX(*), IPIV(*), X(*)
C
      K1 = MAX (N,L)
      IF     (IABS (IER).EQ.N)         GO TO     5
C--
C--      COMPLETE HOUSEHOLDER TRANSFORMATION IF IER LESS THAN N.
C
      IST = IER * (N+1) + 1
      KS  = IER + 1
C
      DO           4         K=KS,N
      LV = M - K + 1
C
C--      GENERATE VECTOR UK AND TRANSFORMATION PARAMETER BETA.
C
      CALL TLUK (A(IST),N,LV,SIG,BETA)
      J = K1 + K
      AUX(J) = -SIG
      IF     (K.EQ.N)        GO TO     4
C
C--      TRANSFORMATION OF MATRIX A.
C--
      NK = N - K
      IF (LV.EQ.1)                     GO TO     2
      CALL TLSTEP (A(IST),A(IST+1),N,N,LV,NK,BETA)
      GO TO        4
    2 DO           3         J=1,NK
      JST = IST + J
    3 A(JST) = A(JST)*(1.-BETA*A(IST)**2)
      IST = IST + N + 1
    4 IPIV(K) = K
C
C
C--      COMPUTE X FROM A AND DIAGONAL ELEMENTS OF A-TRANSPOSED.
C
    5 DO           40        K=1,N
      KN = N-K+1
      IA = (KN-1)*( N+1)+1
      IK = KN *N
      IL = N*N - K + 1
      II = KN + K1
      PIV=1./AUX(II)
      ID = IPIV(KN)-KN
      JA = IA+1
      JK = IK
      JL = IL
C
      DO           20        J=1,K
      H=0.
      IF     (J.EQ.K .AND. J.LE.N-M1)            H = PIV
      IF     (K.EQ.1)                            GO TO     15
      II = JK
C
      DO           10        I=JA,IK
      II=II+N
   10 H = H-A(I)*X(II)
C
   15 H = H*PIV
      X(JL) = H
      JK = JK - 1
   20 JL = JL - N
C
C--      COMPLETE SYMMETRIC PART.
C
      IF     (K.EQ.1)                  GO TO     40
      JL = IA
      DO           25        J=JA,IK
      JL = JL + N
   25 X(J) = X(JL)
C
C--      INTERCHANGE OF ROWS AND COLUMNS ALREADY FINISHED.
C
      IF     (ID.EQ.0)                 GO TO     40
      DO           30        J=IA,IL,N
      II = J + ID
      H  = X(II)
      X(II) = X(J)
   30 X(J)  = H
C
      ID = ID*N
      DO           35        J=IA,IK
      II = J + ID
      H  = X(II)
      X(II) = X(J)
   35 X(J)  = H
C
   40 CONTINUE
C
      RETURN
      END
