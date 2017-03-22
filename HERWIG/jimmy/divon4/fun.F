*CMZ :          04/09/93  16.19.09  by  Jonathan Butterworth
*-- Author :
      REAL FUNCTION FUN(N,RX)
      DOUBLE PRECISION XDB(15), DFUN, FUNX
      REAL RX(N),SNGL
      EXTERNAL DFUN
      SAVE
      DO 10 I=1,N
      XDB(I)=DBLE(RX(I))
 10   CONTINUE
      FUNX=DFUN(N,XDB)
** JMB
      IF (ABS(FUNX).LT.1.E-36) FUNX=0.D0
      FUN=SNGL(FUNX)
      RETURN
      END
