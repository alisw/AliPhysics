*CMZ :          17/07/98  15.49.05  by  Federico Carminati
*-- Author :
      SUBROUTINE SHDCYS
c	=================

c	Particle decays

*KEEP,LUJETS.
      COMMON /LUJETS/ N,K(200000,5),P(200000,5),V(200000,5)
      SAVE /LUJETS/
*KEEP,SHWATE.
      COMMON /SHWATE/ WEI(200000)

*KEND.
      NOLD = N
      CALL LUEXEC
      DO 100 J = NOLD+1,N
        WEI(J) = WEI(K(J,3))
100	CONTINUE
      RETURN
      END
