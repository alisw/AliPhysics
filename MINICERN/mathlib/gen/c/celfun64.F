*
* $Id$
*
* $Log$
* Revision 1.2  1997/12/15 16:18:35  mclareni
* Changes for the Portland Group f77 compiler inside cpp define CERNLIB_QFPGF77
*
* Revision 1.1.1.1  1996/04/01 15:02:01  mclareni
* Mathlib gen
*
*
#include "gen/pilot.h"
#if defined(CERNLIB_DOUBLE)
      SUBROUTINE WELFUN(W,AK2,SN,CN,DN)
C
#include "gen/imp64.inc"
C
      CHARACTER*(*) NAME
      PARAMETER(NAME='CELFUN/WELFUN')
#endif
#if !defined(CERNLIB_DOUBLE)
      SUBROUTINE CELFUN(W,AK2,SN,CN,DN)
C
      CHARACTER*(*) NAME
      PARAMETER(NAME='CELFUN')
#endif
C
C     Jacobian Elliptic Functions SN(W,M),CN(W,M),DN(W,M) for
C     complex argument  w = u+i*v.
C     Iterates for parameter m = k**2 or mc= 1-m if mc < m to obtain
C     fastest convergence and uses finally Jacobi's
C     imaginary transformation to go back to sn, cn, dn of m.
C
#include "gen/defc64.inc"
     +  W,SN,CN,DN,I
      DIMENSION C(4)

C     MACHINE-DEPENDENT: EPS1=2**-(MB/2), EPS2=2**-(MB+3)
C     Where M = Number of bits in mantissa

      PARAMETER (MB = 64)
      PARAMETER (Z0 = 0, Z1 = 1, Z2 = 2, HF = Z1/2, QU = Z1/4)
#if !defined(CERNLIB_CRAY)
      PARAMETER (I = (0D0,1D0))
#endif
#if defined(CERNLIB_CRAY)
      PARAMETER (I = (0E0,1E0))
#endif
      PARAMETER (PI = 3.14159 26535 89793 24D0, PIH = PI/2)
      PARAMETER (EPS1 = Z2**(-MB/2), EPS2 = Z2**(-(MB+3)))

      CHARACTER*80 ERRTXT

#if defined(CERNLIB_QFPGF77)
      DATA AM0 /-1D20/

      SAVE AM0,C,A,L,BIGK
#endif

#if defined(CERNLIB_QF2C)
#include "gen/gcmpfun.inc"
#endif

#if ! defined(CERNLIB_QFPGF77)
      DATA AM0 /-1D20/

      SAVE AM0,C,A,L,BIGK
#endif

#if !defined(CERNLIB_QF2C)
#include "gen/gcmpfun.inc"
#endif

      AM=ABS(AK2)
      IF(AM .GT. 1) THEN
       WRITE(ERRTXT,101) AM
       CALL MTLPRT(NAME,'C320.1',ERRTXT)
       RETURN
      ENDIF
      IF(AM .LE. HF) THEN
       U=W
       V=GIMAG(W)
       IF(AM .EQ. AM0) GO TO 1
       XK2=AM
       B=SQRT(1-AM)
      ELSE
       U=GIMAG(W)
       V=-W
       IF(AM .EQ. AM0) GO TO 1
       XK2=1-AM
       B=SQRT(AM)
      ENDIF

      AM0=AM
      A=1
      L=4
      C(L)=QU*XK2

C     Gaussian arithmetic-geometric mean. Skipped if previous M.

    2 IF(C(L) .GE. EPS1) THEN
       L=L-1
       C(L)=(QU*(A-B))**2
       A1=HF*(A+B)
       B=SQRT(A*B)
       A=A1
       GO TO 2
      ENDIF
      A=HF*(A+B)
      BIGK=PIH/A

C     Descending Landen-Gauss Trafo for real U

    1 X=SIN(A*U)
      IF(V .NE. 0 .OR. X .NE. 0) THEN
       IF(X .NE. 0) THEN
        X=A/X
        DO 3 J = L,4
        X1=C(J)/X
   3    X=X1+X
        SU=1/X
        DU=1-2*X1*SU
        CU=SIGN(SQRT(ABS(1-SU**2)),BIGK-ABS(U))
       ENDIF
       IF(V .NE. 0) THEN
        Y=A/SINH(A*V)
        DO 4 J = L,4
        Y1=C(J)/Y
        Y=Y-Y1
        IF(Y .EQ. 0) Y=EPS2
    4   CONTINUE
        SV=1/Y
        Y1=2*Y1*SV
        DV=1+Y1
        CV=SIGN(SQRT((SV+Y)*SV),Y1)

C     Evaluation of complex sn, cn, dn from real and imaginary ones and
C     Jacobi's imaginary argument trafo if necessary

        IF(U .NE. 0) THEN
         SS=-SU*SV
         D=1/(XK2*SS**2+1)
         DV=DV*D
         CU=CU*D
         CC=CU*CV
         DD=DU*DV
         DN=GCMPLX(DD,XK2*SS*CC)
         CN=GCMPLX(CC,SS*DD)
         SN=GCMPLX(SU*CV*DV,SV*CU*DU)
        ELSE
         SN=GCMPLX(Z0,SV)
         CN=CV
         DN=DV
        ENDIF
       ELSE
        SN=SU
        CN=CU
        DN=DU
       ENDIF
       IF(AM .GT. HF) THEN
        CN=1/CN
        DN=DN*CN
        SN=I*SN*CN
       ENDIF
      ELSE
       SN=0
       CN=1
       DN=1
      ENDIF
      RETURN
  101 FORMAT('MODULUS AK2 = ',1P,E15.6,' OUT OF RANGE')
      END
