 
C***********************************************************************
 
      SUBROUTINE PYWAUX(IAUX,EPS,WRE,WIM)
 
C...Calculates real and imaginary parts of the auxiliary functions W1
C...and W2; see R. K. Ellis, I. Hinchliffe, M. Soldate and J. J. van
C...der Bij, Nucl. Phys. B297 (1988) 221.
      COMMON/LUDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
      SAVE /LUDAT1/
 
      ASINH(X)=LOG(X+SQRT(X**2+1.))
      ACOSH(X)=LOG(X+SQRT(X**2-1.))
 
      IF(EPS.LT.0.) THEN
        IF(IAUX.EQ.1) WRE=2.*SQRT(1.-EPS)*ASINH(SQRT(-1./EPS))
        IF(IAUX.EQ.2) WRE=4.*(ASINH(SQRT(-1./EPS)))**2
        WIM=0.
      ELSEIF(EPS.LT.1.) THEN
        IF(IAUX.EQ.1) WRE=2.*SQRT(1.-EPS)*ACOSH(SQRT(1./EPS))
        IF(IAUX.EQ.2) WRE=4.*(ACOSH(SQRT(1./EPS)))**2-PARU(1)**2
        IF(IAUX.EQ.1) WIM=-PARU(1)*SQRT(1.-EPS)
        IF(IAUX.EQ.2) WIM=-4.*PARU(1)*ACOSH(SQRT(1./EPS))
      ELSE
        IF(IAUX.EQ.1) WRE=2.*SQRT(EPS-1.)*ASIN(SQRT(1./EPS))
        IF(IAUX.EQ.2) WRE=-4.*(ASIN(SQRT(1./EPS)))**2
        WIM=0.
      ENDIF
 
      RETURN
      END
