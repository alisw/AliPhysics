*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/02/15 17:48:48  mclareni
* Kernlib
*
*
#include "kernnum/pilot.h"
      SUBROUTINE RINV(N,A,IDIM,R,IFAIL)
      REAL R(N),T1,T2,T3
      REAL A(IDIM,N),DET,TEMP,S,
     $     C11,C12,C13,C21,C22,C23,C31,C32,C33
      CHARACTER*6 NAME
      DATA NAME/'RINV'/,KPRNT/0/
C
C     ******************************************************************
C
C     REPLACES A BY ITS INVERSE.
C
C     (PARAMETERS AS FOR REQINV.)
C
C     CALLS ... RFACT, RFINV, F010PR, ABEND.
C
C     ******************************************************************
C
C  TEST FOR PARAMETER ERRORS.
C
      IF((N.LT.1).OR.(N.GT.IDIM)) GO TO 7
C
C  TEST FOR N.LE.3.
C
      IF(N.GT.3) GO TO 6
      IFAIL=0
      IF(N.LT.3) GO TO 4
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
      T1=ABS(A(1,1))
      T2=ABS(A(2,1))
      T3=ABS(A(3,1))
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
    3 IF(DET.EQ.0.) GO TO 8
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
      RETURN
C
    4 IF(N.LT.2) GO TO 5
C
C  N=2 CASE BY CRAMERS RULE.
C
      DET=A(1,1)*A(2,2)-A(1,2)*A(2,1)
      IF(DET.EQ.0.) GO TO 8
      S=1./DET
      C11   =S*A(2,2)
      A(1,2)=-S*A(1,2)
      A(2,1)=-S*A(2,1)
      A(2,2)=S*A(1,1)
      A(1,1)=C11
      RETURN
C
C  N=1 CASE.
C
    5 IF(A(1,1).EQ.0.) GO TO 8
      A(1,1)=1./A(1,1)
      RETURN
C
C  N.GT.3 CASES.  FACTORIZE MATRIX AND INVERT.
C
    6 CALL RFACT(N,A,IDIM,R,IFAIL,DET,JFAIL)
      IF(IFAIL.NE.0) RETURN
      CALL RFINV(N,A,IDIM,R)
      RETURN
C
C  ERROR EXITS.
C
    7 IFAIL=+1
      CALL F010PR(NAME,N,IDIM,K,KPRNT)
      RETURN
C
    8 IFAIL=-1
      RETURN
C
      END
