 
C***********************************************************************
 
      FUNCTION PYSPEN(XREIN,XIMIN,IREIM)
 
C...Calculates real and imaginary part of Spence function; see
C...G. 't Hooft and M. Veltman, Nucl. Phys. B153 (1979) 365.
      COMMON/LUDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
      SAVE /LUDAT1/
      DIMENSION B(0:14)
 
      DATA B/
     & 1.000000E+00,        -5.000000E-01,         1.666667E-01,
     & 0.000000E+00,        -3.333333E-02,         0.000000E+00,
     & 2.380952E-02,         0.000000E+00,        -3.333333E-02,
     & 0.000000E+00,         7.575757E-02,         0.000000E+00,
     &-2.531135E-01,         0.000000E+00,         1.166667E+00/
 
      XRE=XREIN
      XIM=XIMIN
      IF(ABS(1.-XRE).LT.1.E-6.AND.ABS(XIM).LT.1.E-6) THEN
        IF(IREIM.EQ.1) PYSPEN=PARU(1)**2/6.
        IF(IREIM.EQ.2) PYSPEN=0.
        RETURN
      ENDIF
 
      XMOD=SQRT(XRE**2+XIM**2)
      IF(XMOD.LT.1.E-6) THEN
        IF(IREIM.EQ.1) PYSPEN=0.
        IF(IREIM.EQ.2) PYSPEN=0.
        RETURN
      ENDIF
 
      XARG=SIGN(ACOS(XRE/XMOD),XIM)
      SP0RE=0.
      SP0IM=0.
      SGN=1.
      IF(XMOD.GT.1.) THEN
        ALGXRE=LOG(XMOD)
        ALGXIM=XARG-SIGN(PARU(1),XARG)
        SP0RE=-PARU(1)**2/6.-(ALGXRE**2-ALGXIM**2)/2.
        SP0IM=-ALGXRE*ALGXIM
        SGN=-1.
        XMOD=1./XMOD
        XARG=-XARG
        XRE=XMOD*COS(XARG)
        XIM=XMOD*SIN(XARG)
      ENDIF
      IF(XRE.GT.0.5) THEN
        ALGXRE=LOG(XMOD)
        ALGXIM=XARG
        XRE=1.-XRE
        XIM=-XIM
        XMOD=SQRT(XRE**2+XIM**2)
        XARG=SIGN(ACOS(XRE/XMOD),XIM)
        ALGYRE=LOG(XMOD)
        ALGYIM=XARG
        SP0RE=SP0RE+SGN*(PARU(1)**2/6.-(ALGXRE*ALGYRE-ALGXIM*ALGYIM))
        SP0IM=SP0IM-SGN*(ALGXRE*ALGYIM+ALGXIM*ALGYRE)
        SGN=-SGN
      ENDIF
 
      XRE=1.-XRE
      XIM=-XIM
      XMOD=SQRT(XRE**2+XIM**2)
      XARG=SIGN(ACOS(XRE/XMOD),XIM)
      ZRE=-LOG(XMOD)
      ZIM=-XARG
 
      SPRE=0.
      SPIM=0.
      SAVERE=1.
      SAVEIM=0.
      DO 100 I=0,14
      IF(MAX(ABS(SAVERE),ABS(SAVEIM)).LT.1E-30) GOTO 110
      TERMRE=(SAVERE*ZRE-SAVEIM*ZIM)/FLOAT(I+1)
      TERMIM=(SAVERE*ZIM+SAVEIM*ZRE)/FLOAT(I+1)
      SAVERE=TERMRE
      SAVEIM=TERMIM
      SPRE=SPRE+B(I)*TERMRE
      SPIM=SPIM+B(I)*TERMIM
  100 CONTINUE
 
  110 IF(IREIM.EQ.1) PYSPEN=SP0RE+SGN*SPRE
      IF(IREIM.EQ.2) PYSPEN=SP0IM+SGN*SPIM
 
      RETURN
      END
