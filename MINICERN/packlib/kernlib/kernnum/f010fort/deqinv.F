*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/02/15 17:48:48  mclareni
* Kernlib
*
*
#include "kernnum/pilot.h"
      SUBROUTINE DEQINV(N,A,IDIM,R,IFAIL,K,B)
      REAL R(N),T1,T2,T3
      DOUBLE PRECISION A(IDIM,N),B(IDIM,K),DET,TEMP,S,
     $                 B1,B2,C11,C12,C13,C21,C22,C23,C31,C32,C33
      CHARACTER*6 NAME
      DATA NAME/'DEQINV'/,KPRNT/1/
C
C     ******************************************************************
C
C     REPLACES B BY THE SOLUTION X OF A*X=B, AND REPLACES A BY ITS IN-
C     VERSE.
C
C     N            ORDER OF THE SQUARE MATRIX IN ARRAY A.
C
C     A            (DOUBLE PRECISION) TWO-DIMENSIONAL ARRAY CONTAINING
C                  AN N BY N MATRIX.
C
C     IDIM         FIRST DIMENSION PARAMETER OF ARRAYS A AND B.
C
C     R            (REAL) WORKING VECTOR OF LENGTH NOT LESS THAN N.
C
C     IFAIL        OUTPUT PARAMETER.   IFAIL= 0 ... NORMAL EXIT.
C                                      IFAIL=-1 ... SINGULAR MATRIX.
C
C     K            NUMBER OF COLUMNS OF THE MATRIX IN ARRAY B.
C
C     B            (DOUBLE PRECISION) TWO-DIMENSIONAL ARRAY CONTAINING
C                  AN N BY K MATRIX.
C
C     CALLS ... DFACT, DFINV, F010PR, ABEND.
C
C     ******************************************************************
C
C  TEST FOR PARAMETER ERRORS.
C
      IF((N.LT.1).OR.(N.GT.IDIM).OR.(K.LT.1)) GO TO 10
C
C  TEST FOR N.LE.3.
C
      IF(N.GT.3) GO TO 9
      IFAIL=0
      IF(N.LT.3) GO TO 5
C
C  N=3 CASE.
C
C     COMPUTE COFACTORS.
      C11=A(2,2)*A(3,3)-A(2,3)*A(3,2)
      C12=A(2,3)*A(3,1)-A(2,1)*A(3,3)
      C13=A(2,1)*A(3,2)-A(2,2)*A(3,1)
      C21=A(3,2)*A(1,3)-A(3,3)*A(1,2)
      C22=A(3,3)*A(1,1)-A(3,1)*A(1,3)
      C23=A(3,1)*A(1,2)-A(3,2)*A(1,1)
      C31=A(1,2)*A(2,3)-A(1,3)*A(2,2)
      C32=A(1,3)*A(2,1)-A(1,1)*A(2,3)
      C33=A(1,1)*A(2,2)-A(1,2)*A(2,1)
      T1=ABS(SNGL(A(1,1)))
      T2=ABS(SNGL(A(2,1)))
      T3=ABS(SNGL(A(3,1)))
C
C     (SET TEMP=PIVOT AND DET=PIVOT*DET.)
      IF(T1.GE.T2) GO TO 1
         IF(T3.GE.T2) GO TO 2
C        (PIVOT IS A21)
            TEMP=A(2,1)
            DET=C13*C32-C12*C33
            GO TO 3
    1 IF(T3.GE.T1) GO TO 2
C     (PIVOT IS A11)
         TEMP=A(1,1)
         DET=C22*C33-C23*C32
         GO TO 3
C     (PIVOT IS A31)
    2    TEMP=A(3,1)
         DET=C23*C12-C22*C13
C
C     SET ELEMENTS OF INVERSE IN A.
    3 IF(DET.EQ.0D0) GO TO 11
      S=TEMP/DET
      A(1,1)=S*C11
      A(1,2)=S*C21
      A(1,3)=S*C31
      A(2,1)=S*C12
      A(2,2)=S*C22
      A(2,3)=S*C32
      A(3,1)=S*C13
      A(3,2)=S*C23
      A(3,3)=S*C33
C
C     REPLACE B BY AINV*B.
      DO 4 J=1,K
         B1=B(1,J)
         B2=B(2,J)
         B(1,J)=A(1,1)*B1+A(1,2)*B2+A(1,3)*B(3,J)
         B(2,J)=A(2,1)*B1+A(2,2)*B2+A(2,3)*B(3,J)
         B(3,J)=A(3,1)*B1+A(3,2)*B2+A(3,3)*B(3,J)
    4 CONTINUE
      RETURN
C
    5 IF(N.LT.2) GO TO 7
C
C  N=2 CASE BY CRAMERS RULE.
C
      DET=A(1,1)*A(2,2)-A(1,2)*A(2,1)
      IF(DET.EQ.0D0) GO TO 11
      S=1D0/DET
      C11   =S*A(2,2)
      A(1,2)=-S*A(1,2)
      A(2,1)=-S*A(2,1)
      A(2,2)=S*A(1,1)
      A(1,1)=C11
      DO 6 J=1,K
         B1=B(1,J)
         B(1,J)=C11*B1+A(1,2)*B(2,J)
         B(2,J)=A(2,1)*B1+A(2,2)*B(2,J)
    6 CONTINUE
      RETURN
C
C  N=1 CASE.
C
    7 IF(A(1,1).EQ.0D0) GO TO 11
      A(1,1)=1D0/A(1,1)
      DO 8 J=1,K
         B(1,J)=A(1,1)*B(1,J)
    8 CONTINUE
      RETURN
C
C  N.GT.3 CASES.  FACTORIZE MATRIX, INVERT AND SOLVE SYSTEM.
C
    9 CALL DFACT(N,A,IDIM,R,IFAIL,DET,JFAIL)
      IF(IFAIL.NE.0) RETURN
      CALL DFEQN(N,A,IDIM,R,K,B)
      CALL DFINV(N,A,IDIM,R)
      RETURN
C
C  ERROR EXITS.
C
   10 IFAIL=+1
      CALL F010PR(NAME,N,IDIM,K,KPRNT)
      RETURN
C
   11 IFAIL=-1
      RETURN
C
      END
