*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/02/15 17:48:49  mclareni
* Kernlib
*
*
#include "kernnum/pilot.h"
      SUBROUTINE CEQN(N,A,IDIM,R,IFAIL,K,B)
      REAL R(N),T1,T2,T3
      COMPLEX A(IDIM,N),B(IDIM,K),ONE,DET,S,TEMP,
     $        B1,Y1,Y2,L11,L21,L22,L31,L32,L33,U12,U13,U23
      CHARACTER*6 NAME
      DATA NAME/'CEQN'/,KPRNT/1/
      DATA ONE/(1.0,0.0)/
C
C     ******************************************************************
C
C     REPLACES B BY THE SOLUTION X OF A*X=B, AFTER WHICH A IS UNDEFINED.
C
C     (PARAMETERS AS FOR CEQINV.)
C
C     CALLS ... CFACT, CFEQN, F010PR, ABEND.
C
C     ******************************************************************
C
C  TEST FOR PARAMETER ERRORS.
C
      IF((N.LT.1).OR.(N.GT.IDIM).OR.(K.LT.1)) GO TO 11
C
C  TEST FOR N.LE.3.
C
      IF(N.GT.3) GO TO 10
      IFAIL=0
      IF(N.LT.3) GO TO 6
C
C  N=3 CASE.
C
C     FACTORIZE MATRIX A=L*U.
C     (FIRST PIVOT SEARCH)
      T1=ABS(REAL(A(1,1)))+ABS(AIMAG(A(1,1)))
      T2=ABS(REAL(A(2,1)))+ABS(AIMAG(A(2,1)))
      T3=ABS(REAL(A(3,1)))+ABS(AIMAG(A(3,1)))
      IF(T1.GE.T2) GO TO 1
         IF(T3.GE.T2) GO TO 2
C        (PIVOT IS A21)
            M1=2
            M2=1
            M3=3
            GO TO 3
    1 IF(T3.GE.T1) GO TO 2
C     (PIVOT IS A11)
         M1=1
         M2=2
         M3=3
         GO TO 3
C     (PIVOT IS A31)
    2    M1=3
         M2=2
         M3=1
    3 TEMP=A(M1,1)
      IF( REAL(TEMP).EQ.0. .AND. AIMAG(TEMP).EQ.0. ) GO TO 10
      L11=ONE/TEMP
      U12=L11*A(M1,2)
      U13=L11*A(M1,3)
      L22=A(M2,2)-A(M2,1)*U12
      L32=A(M3,2)-A(M3,1)*U12
C     (SECOND PIVOT SEARCH)
      T2=ABS(REAL(L22))+ABS(AIMAG(L22))
      T3=ABS(REAL(L32))+ABS(AIMAG(L32))
      IF(T2.GE.T3) GO TO 4
         I=M2
         M2=M3
         M3=I
         TEMP=L22
         L22=L32
         L32=TEMP
    4 L21=A(M2,1)
      L31=A(M3,1)
      IF( REAL(L22).EQ.0. .AND. AIMAG(L22).EQ.0. ) GO TO 10
      L22=ONE/L22
      U23=L22*(A(M2,3)-L21*U13)
      TEMP=A(M3,3)-L31*U13-L32*U23
      IF( REAL(TEMP).EQ.0. .AND. AIMAG(TEMP).EQ.0. ) GO TO 10
      L33=ONE/TEMP
C
C     SOLVE L*Y=B AND U*X=Y.
      DO 5 J=1,K
         Y1=L11*B(M1,J)
         Y2=L22*(B(M2,J)-L21*Y1)
         B(3,J)=L33*(B(M3,J)-L31*Y1-L32*Y2)
         B(2,J)=Y2-U23*B(3,J)
         B(1,J)=Y1-U12*B(2,J)-U13*B(3,J)
    5 CONTINUE
      RETURN
C
    6 IF(N.LT.2) GO TO 8
C
C  N=2 CASE BY CRAMERS RULE.
C
      DET=A(1,1)*A(2,2)-A(1,2)*A(2,1)
      IF( REAL(DET).EQ.0. .AND. AIMAG(DET).EQ.0. ) GO TO 12
      S=ONE/DET
      DO 7 J=1,K
         B1=B(1,J)
         B(1,J)=S*(A(2,2)*B1-A(1,2)*B(2,J))
         B(2,J)=S*(-A(2,1)*B1+A(1,1)*B(2,J))
    7 CONTINUE
      RETURN
C
C  N=1 CASE.
C
    8 IF( REAL(A(1,1)).EQ.0. .AND. AIMAG(A(1,1)).EQ.0. ) GO TO 12
      S=ONE/A(1,1)
      DO 9 J=1,K
         B(1,J)=S*B(1,J)
    9 CONTINUE
      RETURN
C
C  N.GT.3 CASES.  FACTORIZE MATRIX AND SOLVE SYSTEM.
C
   10 CALL CFACT(N,A,IDIM,R,IFAIL,DET,JFAIL)
      IF(IFAIL.NE.0) RETURN
      CALL CFEQN(N,A,IDIM,R,K,B)
      RETURN
C
C  ERROR EXITS.
C
   11 IFAIL=+1
      CALL F010PR(NAME,N,IDIM,K,KPRNT)
      RETURN
C
   12 IFAIL=-1
      RETURN
C
      END
