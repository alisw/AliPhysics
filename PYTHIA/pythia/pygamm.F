 
C*********************************************************************
 
      FUNCTION PYGAMM(X)
 
C...Gives ordinary Gamma function Gamma(x) for positive, real arguments;
C...see M. Abramowitz, I. A. Stegun: Handbook of Mathematical Functions
C...(Dover, 1965) 6.1.36.
      DIMENSION B(8)
      DATA B/-0.577191652,0.988205891,-0.897056937,0.918206857,
     &-0.756704078,0.482199394,-0.193527818,0.035868343/
 
      NX=INT(X)
      DX=X-NX
 
      PYGAMM=1.
      DXP=1.
      DO 100 I=1,8
      DXP=DXP*DX
      PYGAMM=PYGAMM+B(I)*DXP
  100 CONTINUE
      IF(X.LT.1.) THEN
        PYGAMM=PYGAMM/X
      ELSE
        DO 110 IX=1,NX-1
        PYGAMM=(X-IX)*PYGAMM
  110   CONTINUE
      ENDIF
 
      RETURN
      END
