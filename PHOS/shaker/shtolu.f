*CMZ :          17/07/98  15.49.05  by  Federico Carminati
*-- Author :
      SUBROUTINE SHTOLU(KF,Y,PT,PHI,W)
c	================================

c	Load Particle in /LUJETS/

*KEEP,LUJETS.
      COMMON /LUJETS/ N,K(200000,5),P(200000,5),V(200000,5)
      SAVE /LUJETS/
*KEEP,SHWATE.
      COMMON /SHWATE/ WEI(200000)

*KEND.
      N=N+1
      K(N,1) = 1
      K(N,2) = KF
      K(N,3) = 0
      K(N,4) = 0
      K(N,5) = 0

      AM = ULMASS(KF)
      AMT = SQRT(AM**2+PT**2)
      P(N,1) = PT*COS(PHI)
      P(N,2) = PT*SIN(PHI)
      P(N,3) = AMT*SINH(Y)
      P(N,4) = AMT*COSH(Y)
      P(N,5) = AM
      DO 50 JJ=1,5
        V(N,JJ) = 0.
50	CONTINUE
      WEI(N) = W

      RETURN
      END
