 
C********************************************************************* 
 
      FUNCTION RLU(IDUMMY) 
 
C...Purpose: to generate random numbers uniformly distributed between 
C...0 and 1, excluding the endpoints. 
      COMMON/LUDATR/MRLU(6),RRLU(100) 
      SAVE /LUDATR/ 
      EQUIVALENCE (MRLU1,MRLU(1)),(MRLU2,MRLU(2)),(MRLU3,MRLU(3)), 
     &(MRLU4,MRLU(4)),(MRLU5,MRLU(5)),(MRLU6,MRLU(6)), 
     &(RRLU98,RRLU(98)),(RRLU99,RRLU(99)),(RRLU00,RRLU(100)) 
 
C...Initialize generation from given seed. 
      IF(MRLU2.EQ.0) THEN 
        IJ=MOD(MRLU1/30082,31329) 
        KL=MOD(MRLU1,30082) 
        I=MOD(IJ/177,177)+2 
        J=MOD(IJ,177)+2 
        K=MOD(KL/169,178)+1 
        L=MOD(KL,169) 
        DO 110 II=1,97 
        S=0. 
        T=0.5 
        DO 100 JJ=1,24 
        M=MOD(MOD(I*J,179)*K,179) 
        I=J 
        J=K 
        K=M 
        L=MOD(53*L+1,169) 
        IF(MOD(L*M,64).GE.32) S=S+T 
        T=0.5*T 
  100   CONTINUE 
        RRLU(II)=S 
  110   CONTINUE 
        TWOM24=1. 
        DO 120 I24=1,24 
        TWOM24=0.5*TWOM24 
  120   CONTINUE 
        RRLU98=362436.*TWOM24 
        RRLU99=7654321.*TWOM24 
        RRLU00=16777213.*TWOM24 
        MRLU2=1 
        MRLU3=0 
        MRLU4=97 
        MRLU5=33 
      ENDIF 
 
C...Generate next random number. 
  130 RUNI=RRLU(MRLU4)-RRLU(MRLU5) 
      IF(RUNI.LT.0.) RUNI=RUNI+1. 
      RRLU(MRLU4)=RUNI 
      MRLU4=MRLU4-1 
      IF(MRLU4.EQ.0) MRLU4=97 
      MRLU5=MRLU5-1 
      IF(MRLU5.EQ.0) MRLU5=97 
      RRLU98=RRLU98-RRLU99 
      IF(RRLU98.LT.0.) RRLU98=RRLU98+RRLU00 
      RUNI=RUNI-RRLU98 
      IF(RUNI.LT.0.) RUNI=RUNI+1. 
      IF(RUNI.LE.0.OR.RUNI.GE.1.) GOTO 130 
 
C...Update counters. Random number to output. 
      MRLU3=MRLU3+1 
      IF(MRLU3.EQ.1000000000) THEN 
        MRLU2=MRLU2+1 
        MRLU3=0 
      ENDIF 
      RLU=RUNI 
 
      RETURN 
      END 
