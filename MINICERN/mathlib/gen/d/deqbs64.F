*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/04/01 15:02:17  mclareni
* Mathlib gen
*
*
#include "gen/pilot.h"
#if !defined(CERNLIB_DOUBLE)
      SUBROUTINE RDEQBS(N,XA,XZ,Y,H0,EPS,SUB,W)
#endif
#if defined(CERNLIB_DOUBLE)
      SUBROUTINE DDEQBS(N,XA,XZ,Y,H0,EPS,SUB,W)
#include "gen/imp64.inc"
#endif

C     This subroutine is based on the Algol procedure  diffsys  as
C     presented in R. Bulirsch and J. Stoer, Numerical Treatment of
C     Ordinary Differential Equations by Extrapolation Methods,
C     Numer. Math. 8 (1966) 1-13.  The adaption for integration over
C     a given interval (not only over one step) is due to G. Janin.

      CHARACTER NAME*(*)
      CHARACTER*80 ERRTXT
#if !defined(CERNLIB_DOUBLE)
      PARAMETER (NAME = 'RDEQBS')
#endif
#if defined(CERNLIB_DOUBLE)
      PARAMETER (NAME = 'DDEQBS')
#endif
      LOGICAL LCV,LBO,LBH,LFN

      DIMENSION Y(*),W(N,*)

      PARAMETER (DELTA = 1D-14)
      PARAMETER (Z1 = 1, HF = Z1/2, C1 = 3*Z1/2)
      PARAMETER (C9 = 9, C6 = C9/15)
      PARAMETER (C2 = 16/C9, C3 = 64/C9, C4 = 256/C9, C5 = C9/4)

#if !defined(CERNLIB_DOUBLE)
      ENTRY DEQBS(N,XA,XZ,Y,H0,EPS,SUB,W)
#endif

      IF(N .LT. 1 .OR. XA .EQ. XZ .OR. H0 .EQ. 0) RETURN
      DELTAX=DELTA*ABS(XZ-XA)
      X=XA
      H1=SIGN(ABS(H0),XZ-XA)
      SGH=SIGN(Z1,H1)

   12 DO 1 I = 1,N
      W(I,28)=0
      W(I,36)=0
      W(I,23)=Y(I)
      DO 1 K = 1,6
    1 W(I,K)=0
      IF(SGH*(X+H1-XZ) .LT. 0) THEN
       HH=H1
       LFN=.FALSE.
      ELSE
       HH=XZ-X
       IF(ABS(HH) .LT. DELTAX) RETURN
       LFN=.TRUE.
      END IF
      CALL SUB(X,Y,W(1,27))
      LBH=.FALSE.

    2 IF(ABS(HH) .LT. DELTAX) THEN
       WRITE(ERRTXT,101) X
       CALL MTLPRT(NAME,'D201.1',ERRTXT)
       RETURN
      END IF
      A=X+HH
      FC=C1
      LBO=.FALSE.
      M=1
      IR=2
      IS=3
      JJ=6

      DO 3 J = 0,9
      IF(LBO) THEN
       W(1,30)=C2
       W(1,32)=C3
       W(1,34)=C4
      ELSE
       W(1,30)=C5
       W(1,32)=C9
       W(1,34)=36
      END IF
      LCV=J .GT. 2
      IF(J .GT. 6) THEN
       L=6
       W(1,35)=64
       FC=C6*FC
      ELSE
       L=J
       W(1,L+29)=M*M
      END IF
      M=M+M
      G=HH/M
      B=G+G
      IF(LBH .AND. J .LT. 8) THEN
       DO 4 I = 1,N
       W(I,25)=W(I,J+15)
    4  W(I,24)=W(I,J+7)
      ELSE
       KK=(M-2)/2
       M=M-1
       DO 5 I = 1,N
       W(I,24)=W(I,23)
    5  W(I,25)=W(I,23)+G*W(I,27)
       DO 6 K = 1,M
       CALL SUB(X+K*G,W(1,25),W(1,26))
       DO 7 I = 1,N
       U=W(I,24)+B*W(I,26)
       W(I,24)=W(I,25)
       W(I,25)=U
    7  W(I,28)=MAX(W(I,28),ABS(U))
       IF(K .EQ. KK .AND. K .NE. 2) THEN
        JJ=JJ+1
        DO 8 I = 1,N
        W(I,JJ+8)=W(I,25)
    8   W(I,JJ)=W(I,24)
       END IF
    6  CONTINUE
      END IF

      CALL SUB(A,W(1,25),W(1,26))
      DO 9 I = 1,N
      V=W(I,36)
      W(I,36)=HF*(W(I,25)+W(I,24)+G*W(I,26))
      C=W(I,36)
      TA=C
      DO 10 K = 1,L
      B1=W(1,K+29)*V
      B=B1-C
      U=V
      IF(B .NE. 0) THEN
       B=(C-V)/B
       U=C*B
       C=B1*B
      END IF
      V=W(I,K)
      W(I,K)=U
   10 TA=U+TA
      IF(ABS(Y(I)-TA) .GT. EPS*W(I,28)) LCV=.FALSE.
    9 Y(I)=TA
      IF(LCV) THEN
       X=A
       H0=H1
       IF(LFN .OR. ABS(X-XZ) .LT. DELTAX) RETURN
       H1=FC*H1
       GOTO 12
      END IF
      W(1,31)=4
      W(1,33)=16
      LBO=.NOT.LBO
      M=IR
      IR=IS
      IS=M+M
    3 CONTINUE

      LBH=.NOT.LBH
      HH=HF*HH
      H1=HH
      LFN=.FALSE.
      GO TO 2
  101 FORMAT('TOO HIGH ACCURACY REQUIRED NEAR  X = ',D15.8)
      END
