 
C*********************************************************************
 
      SUBROUTINE PYGVMD(ISET,KF,X,Q2,P2,ALAM,XPGA)
C...Purpose: to evaluate the VMD parton distributions of a photon,
C...evolved homogeneously from an initial scale P2 to Q2.
C...Does not include dipole suppression factor.
C...ISET is parton distribution set, see above;
C...additionally ISET=0 is used for the evolution of an anomalous photon
C...which branched at a scale P2 and then evolved homogeneously to Q2.
C...ALAM is the 4-flavour Lambda, which is automatically converted
C...to 3- and 5-flavour equivalents as needed.
      DIMENSION XPGA(-6:6)
      DATA PMC/1.3/, PMB/4.6/, AEM/0.007297/, AEM2PI/0.0011614/
 
C...Reset output.
      DO 100 KFL=-6,6
      XPGA(KFL)=0.
  100 CONTINUE
      KFA=IABS(KF)
 
C...Calculate Lambda; protect against unphysical Q2 and P2 input.
      ALAM3=ALAM*(PMC/ALAM)**(2./27.)
      ALAM5=ALAM*(ALAM/PMB)**(2./23.)
      P2EFF=MAX(P2,1.2*ALAM3**2)
      IF(KFA.EQ.4) P2EFF=MAX(P2EFF,PMC**2)
      IF(KFA.EQ.5) P2EFF=MAX(P2EFF,PMB**2)
      Q2EFF=MAX(Q2,P2EFF)
 
C...Find number of flavours at lower and upper scale.
      NFP=4
      IF(P2EFF.LT.PMC**2) NFP=3
      IF(P2EFF.GT.PMB**2) NFP=5
      NFQ=4
      IF(Q2EFF.LT.PMC**2) NFQ=3
      IF(Q2EFF.GT.PMB**2) NFQ=5
 
C...Find s as sum of 3-, 4- and 5-flavour parts.
      S=0.
      IF(NFP.EQ.3) THEN
        Q2DIV=PMC**2
        IF(NFQ.EQ.3) Q2DIV=Q2EFF
        S=S+(6./27.)*LOG(LOG(Q2DIV/ALAM3**2)/LOG(P2EFF/ALAM3**2))
      ENDIF
      IF(NFP.LE.4.AND.NFQ.GE.4) THEN
        P2DIV=P2EFF
        IF(NFP.EQ.3) P2DIV=PMC**2
        Q2DIV=Q2EFF
        IF(NFQ.EQ.5) Q2DIV=PMB**2
        S=S+(6./25.)*LOG(LOG(Q2DIV/ALAM**2)/LOG(P2DIV/ALAM**2))
      ENDIF
      IF(NFQ.EQ.5) THEN
        P2DIV=PMB**2
        IF(NFP.EQ.5) P2DIV=P2EFF
        S=S+(6./23.)*LOG(LOG(Q2EFF/ALAM5**2)/LOG(P2DIV/ALAM5**2))
      ENDIF
 
C...Calculate frequent combinations of x and s.
      X1=1.-X
      XL=-LOG(X)
      S2=S**2
      S3=S**3
      S4=S**4
 
C...Evaluate homogeneous anomalous parton distributions below or
C...above threshold.
      IF(ISET.EQ.0) THEN
      IF(Q2.LE.P2.OR.(KFA.EQ.4.AND.Q2.LT.PMC**2).OR.
     &(KFA.EQ.5.AND.Q2.LT.PMB**2)) THEN
        XVAL = X * 1.5 * (X**2+X1**2)
        XGLU = 0.
        XSEA = 0.
      ELSE
        XVAL = (1.5/(1.-0.197*S+4.33*S2)*X**2 + (1.5+2.10*S)/
     &  (1.+3.29*S)*X1**2 + 5.23*S/(1.+1.17*S+19.9*S3)*X*X1) *
     &  X**(1./(1.+1.5*S)) * (1.-X**2)**(2.667*S)
        XGLU = 4.*S/(1.+4.76*S+15.2*S2+29.3*S4) *
     &  X**(-2.03*S/(1.+2.44*S)) * (X1*XL)**(1.333*S) *
     &  ((4.*X**2+7.*X+4.)*X1/3. - 2.*X*(1.+X)*XL)
        XSEA = S2/(1.+4.54*S+8.19*S2+8.05*S3) *
     &  X**(-1.54*S/(1.+1.29*S)) * X1**(2.667*S) *
     &  ((8.-73.*X+62.*X**2)*X1/9. + (3.-8.*X**2/3.)*X*XL +
     &  (2.*X-1.)*X*XL**2)
      ENDIF
 
C...Evaluate set 1D parton distributions below or above threshold.
      ELSEIF(ISET.EQ.1) THEN
      IF(Q2.LE.P2.OR.(KFA.EQ.4.AND.Q2.LT.PMC**2).OR.
     &(KFA.EQ.5.AND.Q2.LT.PMB**2)) THEN
        XVAL = 1.294 * X**0.80 * X1**0.76
        XGLU = 1.273 * X**0.40 * X1**1.76
        XSEA = 0.100 * X1**3.76
      ELSE
        XVAL = 1.294/(1.+0.252*S+3.079*S2) * X**(0.80-0.13*S) *
     &  X1**(0.76+0.667*S) * XL**(2.*S)
        XGLU = 7.90*S/(1.+5.50*S) * EXP(-5.16*S) *
     &  X**(-1.90*S/(1.+3.60*S)) * X1**1.30 * XL**(0.50+3.*S) +
     &  1.273 * EXP(-10.*S) * X**0.40 * X1**(1.76+3.*S)
        XSEA = (0.1-0.397*S2+1.121*S3)/(1.+5.61*S2+5.26*S3) *
     &  X**(-7.32*S2/(1.+10.3*S2)) *
     &  X1**((3.76+15.*S+12.*S2)/(1.+4.*S))
        XSEA0 = 0.100 * X1**3.76
      ENDIF
 
C...Evaluate set 1M parton distributions below or above threshold.
      ELSEIF(ISET.EQ.2) THEN
      IF(Q2.LE.P2.OR.(KFA.EQ.4.AND.Q2.LT.PMC**2).OR.
     &(KFA.EQ.5.AND.Q2.LT.PMB**2)) THEN
        XVAL = 0.8477 * X**0.51 * X1**1.37
        XGLU = 3.42 * X**0.255 * X1**2.37
        XSEA = 0.
      ELSE
        XVAL = 0.8477/(1.+1.37*S+2.18*S2+3.73*S3) * X**(0.51+0.21*S)
     &  * X1**1.37 * XL**(2.667*S)
        XGLU = 24.*S/(1.+9.6*S+0.92*S2+14.34*S3) * EXP(-5.94*S) *
     &  X**((-0.013-1.80*S)/(1.+3.14*S)) * X1**(2.37+0.4*S) *
     &  XL**(0.32+3.6*S) + 3.42 * EXP(-12.*S) * X**0.255 *
     &  X1**(2.37+3.*S)
        XSEA = 0.842*S/(1.+21.3*S-33.2*S2+229.*S3) *
     &  X**((0.13-2.90*S)/(1.+5.44*S)) * X1**(3.45+0.5*S) *
     &  XL**(2.8*S)
        XSEA0 = 0.
      ENDIF
 
C...Evaluate set 2D parton distributions below or above threshold.
      ELSEIF(ISET.EQ.3) THEN
      IF(Q2.LE.P2.OR.(KFA.EQ.4.AND.Q2.LT.PMC**2).OR.
     &(KFA.EQ.5.AND.Q2.LT.PMB**2)) THEN
        XVAL = X**0.46 * X1**0.64 + 0.76 * X
        XGLU = 1.925 * X1**2
        XSEA = 0.242 * X1**4
      ELSE
        XVAL = (1.+0.186*S)/(1.-0.209*S+1.495*S2) * X**(0.46+0.25*S)
     &  * X1**((0.64+0.14*S+5.*S2)/(1.+S)) * XL**(1.9*S) +
     &  (0.76+0.4*S) * X * X1**(2.667*S)
        XGLU = (1.925+5.55*S+147.*S2)/(1.-3.59*S+3.32*S2) *
     &  EXP(-18.67*S) * X**((-5.81*S-5.34*S2)/(1.+29.*S-4.26*S2))
     &  * X1**((2.-5.9*S)/(1.+1.7*S)) * XL**(9.3*S/(1.+1.7*S))
        XSEA = (0.242-0.252*S+1.19*S2)/(1.-0.607*S+21.95*S2) *
     &  X**(-12.1*S2/(1.+2.62*S+16.7*S2)) * X1**4 * XL**S
        XSEA0 = 0.242 * X1**4
      ENDIF
 
C...Evaluate set 2M parton distributions below or above threshold.
      ELSEIF(ISET.EQ.4) THEN
      IF(Q2.LE.P2.OR.(KFA.EQ.4.AND.Q2.LT.PMC**2).OR.
     &(KFA.EQ.5.AND.Q2.LT.PMB**2)) THEN
        XVAL = 1.168 * X**0.50 * X1**2.60 + 0.965 * X
        XGLU = 1.808 * X1**2
        XSEA = 0.209 * X1**4
      ELSE
        XVAL = (1.168+1.771*S+29.35*S2) * EXP(-5.776*S) *
     &  X**((0.5+0.208*S)/(1.-0.794*S+1.516*S2)) *
     &  X1**((2.6+7.6*S)/(1.+5.*S)) * XL**(5.15*S/(1.+2.*S)) +
     &  (0.965+22.35*S)/(1.+18.4*S) * X * X1**(2.667*S)
        XGLU = (1.808+29.9*S)/(1.+26.4*S) * EXP(-5.28*S) *
     &  X**((-5.35*S-10.11*S2)/(1.+31.71*S)) *
     &  X1**((2.-7.3*S+4.*S2)/(1.+2.5*S)) *
     &  XL**(10.9*S/(1.+2.5*S))
        XSEA = (0.209+0.644*S2)/(1.+0.319*S+17.6*S2) *
     &  X**((-0.373*S-7.71*S2)/(1.+0.815*S+11.0*S2)) *
     &  X1**(4.+S) * XL**(0.45*S)
        XSEA0 = 0.209 * X1**4
      ENDIF
      ENDIF
 
C...Threshold factors for c and b sea.
      SLL=LOG(LOG(Q2EFF/ALAM**2)/LOG(P2EFF/ALAM**2))
      XCHM=0.
      IF(Q2.GT.PMC**2.AND.Q2.GT.1.001*P2EFF) THEN
        SCH=MAX(0.,LOG(LOG(PMC**2/ALAM**2)/LOG(P2EFF/ALAM**2)))
        IF(ISET.EQ.0) THEN
          XCHM=XSEA*(1.-(SCH/SLL)**2)
        ELSE
          XCHM=MAX(0.,XSEA-XSEA0*X1**(2.667*S))*(1.-SCH/SLL)
        ENDIF
      ENDIF
      XBOT=0.
      IF(Q2.GT.PMB**2.AND.Q2.GT.1.001*P2EFF) THEN
        SBT=MAX(0.,LOG(LOG(PMB**2/ALAM**2)/LOG(P2EFF/ALAM**2)))
        IF(ISET.EQ.0) THEN
          XBOT=XSEA*(1.-(SBT/SLL)**2)
        ELSE
          XBOT=MAX(0.,XSEA-XSEA0*X1**(2.667*S))*(1.-SBT/SLL)
        ENDIF
      ENDIF
 
C...Fill parton distributions.
      XPGA(0)=XGLU
      XPGA(1)=XSEA
      XPGA(2)=XSEA
      XPGA(3)=XSEA
      XPGA(4)=XCHM
      XPGA(5)=XBOT
      XPGA(KFA)=XPGA(KFA)+XVAL
      DO 110 KFL=1,5
      XPGA(-KFL)=XPGA(KFL)
  110 CONTINUE
 
      RETURN
      END
