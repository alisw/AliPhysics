 
C*********************************************************************
 
      SUBROUTINE PYGANO(KF,X,Q2,P2,ALAM,XPGA)
C...Purpose: to evaluate the parton distributions of the anomalous
C...photon, inhomogeneously evolved from a scale P2 (where it vanishes)
C...to Q2.
C...KF=0 gives the sum over (up to) 5 flavours,
C...KF<0 limits to flavours up to abs(KF),
C...KF>0 is for flavour KF only.
C...ALAM is the 4-flavour Lambda, which is automatically converted
C...to 3- and 5-flavour equivalents as needed.
      DIMENSION XPGA(-6:6),ALAMSQ(3:5)
      DATA PMC/1.3/, PMB/4.6/, AEM/0.007297/, AEM2PI/0.0011614/
 
C...Reset output.
      DO 100 KFL=-6,6
      XPGA(KFL)=0.
  100 CONTINUE
      IF(Q2.LE.P2) RETURN
      KFA=IABS(KF)
 
C...Calculate Lambda; protect against unphysical Q2 and P2 input.
      ALAMSQ(3)=(ALAM*(PMC/ALAM)**(2./27.))**2
      ALAMSQ(4)=ALAM**2
      ALAMSQ(5)=(ALAM*(ALAM/PMB)**(2./23.))**2
      P2EFF=MAX(P2,1.2*ALAMSQ(3))
      IF(KF.EQ.4) P2EFF=MAX(P2EFF,PMC**2)
      IF(KF.EQ.5) P2EFF=MAX(P2EFF,PMB**2)
      Q2EFF=MAX(Q2,P2EFF)
      XL=-LOG(X)
 
C...Find number of flavours at lower and upper scale.
      NFP=4
      IF(P2EFF.LT.PMC**2) NFP=3
      IF(P2EFF.GT.PMB**2) NFP=5
      NFQ=4
      IF(Q2EFF.LT.PMC**2) NFQ=3
      IF(Q2EFF.GT.PMB**2) NFQ=5
 
C...Define range of flavour loop.
      IF(KF.EQ.0) THEN
        KFLMN=1
        KFLMX=5
      ELSEIF(KF.LT.0) THEN
        KFLMN=1
        KFLMX=KFA
      ELSE
        KFLMN=KFA
        KFLMX=KFA
      ENDIF
 
C...Loop over flavours the photon can branch into.
      DO 110 KFL=KFLMN,KFLMX
 
C...Light flavours: calculate t range and (approximate) s range.
      IF(KFL.LE.3.AND.(KFL.EQ.1.OR.KFL.EQ.KF)) THEN
        TDIFF=LOG(Q2EFF/P2EFF)
        S=(6./(33.-2.*NFQ))*LOG(LOG(Q2EFF/ALAMSQ(NFQ))/
     &  LOG(P2EFF/ALAMSQ(NFQ)))
        IF(NFQ.GT.NFP) THEN
          Q2DIV=PMB**2
          IF(NFQ.EQ.4) Q2DIV=PMC**2
          SNFQ=(6./(33.-2.*NFQ))*LOG(LOG(Q2DIV/ALAMSQ(NFQ))/
     &    LOG(P2EFF/ALAMSQ(NFQ)))
          SNFP=(6./(33.-2.*(NFQ-1)))*LOG(LOG(Q2DIV/ALAMSQ(NFQ-1))/
     &    LOG(P2EFF/ALAMSQ(NFQ-1)))
          S=S+(LOG(Q2DIV/P2EFF)/LOG(Q2EFF/P2EFF))*(SNFP-SNFQ)
        ENDIF
        IF(NFQ.EQ.5.AND.NFP.EQ.3) THEN
          Q2DIV=PMC**2
          SNF4=(6./(33.-2.*4))*LOG(LOG(Q2DIV/ALAMSQ(4))/
     &    LOG(P2EFF/ALAMSQ(4)))
          SNF3=(6./(33.-2.*3))*LOG(LOG(Q2DIV/ALAMSQ(3))/
     &    LOG(P2EFF/ALAMSQ(3)))
          S=S+(LOG(Q2DIV/P2EFF)/LOG(Q2EFF/P2EFF))*(SNF3-SNF4)
        ENDIF
 
C...u and s quark do not need a separate treatment when d has been done.
      ELSEIF(KFL.EQ.2.OR.KFL.EQ.3) THEN
 
C...Charm: as above, but only include range above c threshold.
      ELSEIF(KFL.EQ.4) THEN
        IF(Q2.LE.PMC**2) GOTO 110
        P2EFF=MAX(P2EFF,PMC**2)
        Q2EFF=MAX(Q2EFF,P2EFF)
        TDIFF=LOG(Q2EFF/P2EFF)
        S=(6./(33.-2.*NFQ))*LOG(LOG(Q2EFF/ALAMSQ(NFQ))/
     &  LOG(P2EFF/ALAMSQ(NFQ)))
        IF(NFQ.EQ.5.AND.NFP.EQ.4) THEN
          Q2DIV=PMB**2
          SNFQ=(6./(33.-2.*NFQ))*LOG(LOG(Q2DIV/ALAMSQ(NFQ))/
     &    LOG(P2EFF/ALAMSQ(NFQ)))
          SNFP=(6./(33.-2.*(NFQ-1)))*LOG(LOG(Q2DIV/ALAMSQ(NFQ-1))/
     &    LOG(P2EFF/ALAMSQ(NFQ-1)))
          S=S+(LOG(Q2DIV/P2EFF)/LOG(Q2EFF/P2EFF))*(SNFP-SNFQ)
        ENDIF
 
C...Bottom: as above, but only include range above b threshold.
      ELSEIF(KFL.EQ.5) THEN
        IF(Q2.LE.PMB**2) GOTO 110
        P2EFF=MAX(P2EFF,PMB**2)
        Q2EFF=MAX(Q2,P2EFF)
        TDIFF=LOG(Q2EFF/P2EFF)
        S=(6./(33.-2.*NFQ))*LOG(LOG(Q2EFF/ALAMSQ(NFQ))/
     &  LOG(P2EFF/ALAMSQ(NFQ)))
      ENDIF
 
C...Evaluate flavour-dependent prefactor (charge^2 etc.).
      CHSQ=1./9.
      IF(KFL.EQ.2.OR.KFL.EQ.4) CHSQ=4./9.
      FAC=AEM2PI*2.*CHSQ*TDIFF
 
C...Evaluate parton distributions (normalized to unit momentum sum).
      IF(KFL.EQ.1.OR.KFL.EQ.4.OR.KFL.EQ.5.OR.KFL.EQ.KF) THEN
        XVAL= ((1.5+2.49*S+26.9*S**2)/(1.+32.3*S**2)*X**2 +
     &  (1.5-0.49*S+7.83*S**2)/(1.+7.68*S**2)*(1.-X)**2 +
     &  1.5*S/(1.-3.2*S+7.*S**2)*X*(1.-X)) *
     &  X**(1./(1.+0.58*S)) * (1.-X**2)**(2.5*S/(1.+10.*S))
        XGLU= 2.*S/(1.+4.*S+7.*S**2) *
     &  X**(-1.67*S/(1.+2.*S)) * (1.-X**2)**(1.2*S) *
     &  ((4.*X**2+7.*X+4.)*(1.-X)/3. - 2.*X*(1.+X)*XL)
        XSEA= 0.333*S**2/(1.+4.90*S+4.69*S**2+21.4*S**3) *
     &  X**(-1.18*S/(1.+1.22*S)) * (1.-X)**(1.2*S) *
     &  ((8.-73.*X+62.*X**2)*(1.-X)/9. + (3.-8.*X**2/3.)*X*XL +
     &  (2.*X-1.)*X*XL**2)
 
C...Threshold factors for c and b sea.
        SLL=LOG(LOG(Q2EFF/ALAM**2)/LOG(P2EFF/ALAM**2))
        XCHM=0.
        IF(Q2.GT.PMC**2.AND.Q2.GT.1.001*P2EFF) THEN
          SCH=MAX(0.,LOG(LOG(PMC**2/ALAM**2)/LOG(P2EFF/ALAM**2)))
          XCHM=XSEA*(1.-(SCH/SLL)**3)
        ENDIF
        XBOT=0.
        IF(Q2.GT.PMB**2.AND.Q2.GT.1.001*P2EFF) THEN
          SBT=MAX(0.,LOG(LOG(PMB**2/ALAM**2)/LOG(P2EFF/ALAM**2)))
          XBOT=XSEA*(1.-(SBT/SLL)**3)
        ENDIF
      ENDIF
 
C...Add contribution of each valence flavour.
      XPGA(0)=XPGA(0)+FAC*XGLU
      XPGA(1)=XPGA(1)+FAC*XSEA
      XPGA(2)=XPGA(2)+FAC*XSEA
      XPGA(3)=XPGA(3)+FAC*XSEA
      XPGA(4)=XPGA(4)+FAC*XCHM
      XPGA(5)=XPGA(5)+FAC*XBOT
      XPGA(KFL)=XPGA(KFL)+FAC*XVAL
  110 CONTINUE
      DO 120 KFL=1,5
      XPGA(-KFL)=XPGA(KFL)
  120 CONTINUE
 
      RETURN
      END
