*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/04/01 15:02:01  mclareni
* Mathlib gen
*
*
#include "gen/pilot.h"
#if defined(CERNLIB_DOUBLE)
      SUBROUTINE DELFUN(U,AK2,SN,CN,DN)
#include "gen/imp64.inc"
#endif
#if defined(CERNLIB_SINGLE)
      SUBROUTINE RELFUN(U,AK2,SN,CN,DN)
#endif
      DIMENSION C(4)

C     Machine-dependent: EPS1=2**-(MB/2), EPS2=2**-(MB+3)
C     Where M = Number of bits in mantissa

      PARAMETER (MB = 64)
      PARAMETER (Z0 = 0, Z1 = 1, Z2 = 2, HF = Z1/2, QU = Z1/4)
      PARAMETER (PI = 3.14159 26535 89793 24D0, PIH = PI/2)
      PARAMETER (EPS1 = Z2**(-MB/2), EPS2 = Z2**(-(MB+3)))

      DATA AM0 /-1D20/

      SAVE AM0,C,A,L,BIGK

#if defined(CERNLIB_SINGLE)
      ENTRY ELFUN(U,AK2,SN,CN,DN)
#endif

      XM=ABS(AK2)
      IF(U .EQ. 0) THEN
       SN=0
       DN=1
       CN=1
      ELSEIF(XM .EQ. 0) THEN
       SN=SIN(U)
       DN=1
       CN=COS(U)
      ELSEIF(XM .EQ. 1) THEN
       SN=TANH(U)
       DN=1/COSH(U)
       CN=DN
      ELSE
       IF(XM .LE. 1) THEN
        U1=U
        AM=XM
       ELSE
        W=SQRT(XM)
        U1=W*U
        AM=1/XM
       ENDIF
       IF(AM .LE. HF) THEN
        IF(AM .EQ. AM0) GO TO 1
        AM0=AM
        C(4)=QU*AM
        B=SQRT(1-AM)
       ELSE
        AMC=1-AM
        IF(AMC .EQ. AM0) GO TO 1
        AM0=AMC
        C(4)=QU*AMC
        B=SQRT(AM)
       ENDIF

C     Gauss arithmetic-geometric mean. Skipped if previous modulus.

       A=1
       L=4
    2  IF(C(L) .GE. EPS1) THEN
        L=L-1
        C(L)=(QU*(A-B))**2
        A1=HF*(A+B)
        B=SQRT(A*B)
        A=A1
        GO TO 2
       ENDIF
       A=HF*(A+B)
       BIGK=PIH/A

C     Descending Landen-Gauss transformation for real argument

    1  IF(AM .LE. HF) THEN
        X=SIN(A*U1)
        IF(X .EQ. 0) THEN
         SN=0
         DN=1
         CN=1
        ELSE
         X=A/X
         DO 3 J = L,4
         X1=C(J)/X
    3    X=X1+X
         H=1/X
         SN=H
         DN=1-2*X1*H
         CN=SIGN(SQRT(ABS(1-H**2)),BIGK-ABS(U1))
        ENDIF
       ELSE

C     Descending Landen-Gauss Transformation for imaginary argument

        Y=A/SINH(A*U1)
        DO 4 J=L,4
        Y1=C(J)/Y
        Y=Y-Y1
    4   IF(Y .EQ. 0) Y=EPS2
        H=1/Y
        Y1=2*Y1*H
        CN=SIGN(SQRT(Y/(H+Y)),Y1)
        DN=CN*(1+Y1)
        SN=CN*H
       ENDIF
       IF(XM .GT. 1) THEN
        SN=SN/W
        H=DN
        DN=CN
        CN=H
       ENDIF
      ENDIF
      RETURN
      END
