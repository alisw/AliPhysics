 
C*********************************************************************
 
      SUBROUTINE PYSTPI(X,Q2,XPPI)
 
C...Gives pi+ structure function according to two different
C...parametrizations.
      COMMON/LUDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
      COMMON/PYPARS/MSTP(200),PARP(200),MSTI(200),PARI(200)
      COMMON/PYINT1/MINT(400),VINT(400)
      SAVE /LUDAT1/
      SAVE /PYPARS/,/PYINT1/
      DIMENSION XPPI(-6:6),COW(3,5,4,2),XQ(9),TS(6)
 
C...The following data lines are coefficients needed in the
C...Owens pion structure function parametrizations, see below.
C...Expansion coefficients for up and down valence quark distributions.
      DATA ((COW(IP,IS,1,1),IS=1,5),IP=1,3)/
     1  4.0000E-01,  7.0000E-01,  0.0000E+00,  0.0000E+00,  0.0000E+00,
     2 -6.2120E-02,  6.4780E-01,  0.0000E+00,  0.0000E+00,  0.0000E+00,
     3 -7.1090E-03,  1.3350E-02,  0.0000E+00,  0.0000E+00,  0.0000E+00/
      DATA ((COW(IP,IS,1,2),IS=1,5),IP=1,3)/
     1  4.0000E-01,  6.2800E-01,  0.0000E+00,  0.0000E+00,  0.0000E+00,
     2 -5.9090E-02,  6.4360E-01,  0.0000E+00,  0.0000E+00,  0.0000E+00,
     3 -6.5240E-03,  1.4510E-02,  0.0000E+00,  0.0000E+00,  0.0000E+00/
C...Expansion coefficients for gluon distribution.
      DATA ((COW(IP,IS,2,1),IS=1,5),IP=1,3)/
     1  8.8800E-01,  0.0000E+00,  3.1100E+00,  6.0000E+00,  0.0000E+00,
     2 -1.8020E+00, -1.5760E+00, -1.3170E-01,  2.8010E+00, -1.7280E+01,
     3  1.8120E+00,  1.2000E+00,  5.0680E-01, -1.2160E+01,  2.0490E+01/
      DATA ((COW(IP,IS,2,2),IS=1,5),IP=1,3)/
     1  7.9400E-01,  0.0000E+00,  2.8900E+00,  6.0000E+00,  0.0000E+00,
     2 -9.1440E-01, -1.2370E+00,  5.9660E-01, -3.6710E+00, -8.1910E+00,
     3  5.9660E-01,  6.5820E-01, -2.5500E-01, -2.3040E+00,  7.7580E+00/
C...Expansion coefficients for (up+down+strange) quark sea distribution.
      DATA ((COW(IP,IS,3,1),IS=1,5),IP=1,3)/
     1  9.0000E-01,  0.0000E+00,  5.0000E+00,  0.0000E+00,  0.0000E+00,
     2 -2.4280E-01, -2.1200E-01,  8.6730E-01,  1.2660E+00,  2.3820E+00,
     3  1.3860E-01,  3.6710E-03,  4.7470E-02, -2.2150E+00,  3.4820E-01/
      DATA ((COW(IP,IS,3,2),IS=1,5),IP=1,3)/
     1  9.0000E-01,  0.0000E+00,  5.0000E+00,  0.0000E+00,  0.0000E+00,
     2 -1.4170E-01, -1.6970E-01, -2.4740E+00, -2.5340E+00,  5.6210E-01,
     3 -1.7400E-01, -9.6230E-02,  1.5750E+00,  1.3780E+00, -2.7010E-01/
C...Expansion coefficients for charm quark sea distribution.
      DATA ((COW(IP,IS,4,1),IS=1,5),IP=1,3)/
     1  0.0000E+00, -2.2120E-02,  2.8940E+00,  0.0000E+00,  0.0000E+00,
     2  7.9280E-02, -3.7850E-01,  9.4330E+00,  5.2480E+00,  8.3880E+00,
     3 -6.1340E-02, -1.0880E-01, -1.0852E+01, -7.1870E+00, -1.1610E+01/
      DATA ((COW(IP,IS,4,2),IS=1,5),IP=1,3)/
     1  0.0000E+00, -8.8200E-02,  1.9240E+00,  0.0000E+00,  0.0000E+00,
     2  6.2290E-02, -2.8920E-01,  2.4240E-01, -4.4630E+00, -8.3670E-01,
     3 -4.0990E-02, -1.0820E-01,  2.0360E+00,  5.2090E+00, -4.8400E-02/
 
C...Euler's beta function, requires ordinary Gamma function
      EULBET(X,Y)=PYGAMM(X)*PYGAMM(Y)/PYGAMM(X+Y)
 
C...Reset output array.
      DO 100 KFL=-6,6
      XPPI(KFL)=0.
  100 CONTINUE
 
      IF(MSTP(53).LE.2) THEN
C...Pion structure functions from Owens.
C...Allowed variable range: 4 GeV^2 < Q^2 < approx 2000 GeV^2.
 
C...Determine set, Lambda and s expansion variable.
        NSET=MSTP(53)
        IF(NSET.EQ.1) ALAM=0.2
        IF(NSET.EQ.2) ALAM=0.4
        VINT(231)=4.
        IF(MSTP(57).LE.0) THEN
          SD=0.
        ELSE
          Q2IN=MIN(2E3,MAX(4.,Q2))
          SD=LOG(LOG(Q2IN/ALAM**2)/LOG(4./ALAM**2))
        ENDIF
 
C...Calculate structure functions.
        DO 120 KFL=1,4
        DO 110 IS=1,5
        TS(IS)=COW(1,IS,KFL,NSET)+COW(2,IS,KFL,NSET)*SD+
     &  COW(3,IS,KFL,NSET)*SD**2
  110   CONTINUE
        IF(KFL.EQ.1) THEN
          XQ(KFL)=X**TS(1)*(1.-X)**TS(2)/EULBET(TS(1),TS(2)+1.)
        ELSE
          XQ(KFL)=TS(1)*X**TS(2)*(1.-X)**TS(3)*(1.+TS(4)*X+TS(5)*X**2)
        ENDIF
  120   CONTINUE
 
C...Put into output array.
        XPPI(0)=XQ(2)
        XPPI(1)=XQ(3)/6.
        XPPI(2)=XQ(1)+XQ(3)/6.
        XPPI(3)=XQ(3)/6.
        XPPI(4)=XQ(4)
        XPPI(-1)=XQ(1)+XQ(3)/6.
        XPPI(-2)=XQ(3)/6.
        XPPI(-3)=XQ(3)/6.
        XPPI(-4)=XQ(4)
 
C...Leading order pion structure functions from Gluck, Reya and Vogt.
C...Allowed variable range: 0.25 GeV^2 < Q^2 < 10^8 GeV^2 and
C...10^-5 < x < 1.
      ELSE
 
C...Determine s expansion variable and some x expressions.
        VINT(231)=0.25
        IF(MSTP(57).LE.0) THEN
          SD=0.
        ELSE
          Q2IN=MIN(1E8,MAX(0.25,Q2))
          SD=LOG(LOG(Q2IN/0.232**2)/LOG(0.25/0.232**2))
        ENDIF
        SD2=SD**2
        XL=-LOG(X)
        XS=SQRT(X)
 
C...Evaluate valence, gluon and sea distributions.
        XFVAL=(0.519+0.180*SD-0.011*SD2)*X**(0.499-0.027*SD)*
     &  (1.+(0.381-0.419*SD)*XS)*(1.-X)**(0.367+0.563*SD)
        XFGLU=(X**(0.482+0.341*SQRT(SD))*((0.678+0.877*SD-0.175*SD2)+
     &  (0.338-1.597*SD)*XS+(-0.233*SD+0.406*SD2)*X)+
     &  SD**0.599*EXP(-(0.618+2.070*SD)+SQRT(3.676*SD**1.263*XL)))*
     &  (1.-X)**(0.390+1.053*SD)
        XFSEA=SD**0.55*(1.-0.748*XS+(0.313+0.935*SD)*X)*(1.-X)**3.359*
     &  EXP(-(4.433+1.301*SD)+SQRT((9.30-0.887*SD)*SD**0.56*XL))/
     &  XL**(2.538-0.763*SD)
        IF(SD.LE.0.888) THEN
          XFCHM=0.
        ELSE
          XFCHM=(SD-0.888)**1.02*(1.+1.008*X)*(1.-X)**(1.208+0.771*SD)*
     &    EXP(-(4.40+1.493*SD)+SQRT((2.032+1.901*SD)*SD**0.39*XL))
        ENDIF
        IF(SD.LE.1.351) THEN
          XFBOT=0.
        ELSE
          XFBOT=(SD-1.351)**1.03*(1.-X)**(0.697+0.855*SD)*
     &    EXP(-(4.51+1.490*SD)+SQRT((3.056+1.694*SD)*SD**0.39*XL))
        ENDIF
 
C...Put into output array.
        XPPI(0)=XFGLU
        XPPI(1)=XFSEA
        XPPI(2)=XFSEA
        XPPI(3)=XFSEA
        XPPI(4)=XFCHM
        XPPI(5)=XFBOT
        DO 130 KFL=1,5
        XPPI(-KFL)=XPPI(KFL)
  130   CONTINUE
        XPPI(2)=XPPI(2)+XFVAL
        XPPI(-1)=XPPI(-1)+XFVAL
      ENDIF
 
      RETURN
      END
