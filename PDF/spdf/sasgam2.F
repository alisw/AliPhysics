#include "pdf/pilot.h"
C...SaSgam version 2 - parton distributions of the photon
C...by Gerhard A. Schuler and Torbjorn Sjostrand
C...For further information see Z. Phys. C68 (1995) 607
C...and CERN-TH/96-04 and LU TP 96-2.
C...Program last changed on 18 January 1996.
 
C!!!Note that one further call parameter - IP2 - has been added
C!!!to the SASGAM argument list compared with version 1.
 
C...The user should only need to call the SASGAM routine,
C...which in turn calls the auxiliary routines SASVMD, SASANO,
C...SASBEH and SASDIR. The package is self-contained.
 
C...One particular aspect of these parametrizations is that F2 for
C...the photon is not obtained just as the charge-squared-weighted
C...sum of quark distributions, but differ in the treatment of
C...heavy flavours (in F2 the DIS relation W2 = Q2*(1-x)/x restricts
C...the kinematics range of heavy-flavour production, but the same
C...kinematics is not relevant e.g. for jet production) and, for the
C...'MSbar' fits, in the addition of a Cgamma term related to the
C...separation of direct processes. Schematically:
C...PDF = VMD (rho, omega, phi) + anomalous (d, u, s, c, b).
C...F2  = VMD (rho, omega, phi) + anomalous (d, u, s) +
C...      Bethe-Heitler (c, b) (+ Cgamma (d, u, s)).
C...The J/psi and Upsilon states have not been included in the VMD sum,
C...but low c and b masses in the other components should compensate
C...for this in a duality sense.
 
C...The calling sequence is the following:
C     CALL SASGAM2(ISET,X,Q2,P2,IP2,F2GM,XPDFGM)
C...with the following declaration statement:
C     DIMENSION XPDFGM(-6:6)
C...and, optionally, further information in:
C     COMMON/SASCOM/XPVMD(-6:6),XPANL(-6:6),XPANH(-6:6),XPBEH(-6:6),
C    &XPDIR(-6:6)
C     COMMON/SASVAL/VXPVMD(-6:6),VXPANL(-6:6),VXPANH(-6:6),VXPDGM(-6:6)
C...Input:  ISET = 1 : SaS set 1D ('DIS',   Q0 = 0.6 GeV)
C                = 2 : SaS set 1M ('MSbar', Q0 = 0.6 GeV)
C                = 3 : SaS set 2D ('DIS',   Q0 =  2  GeV)
C                = 4 : SaS set 2M ('MSbar', Q0 =  2  GeV)
C           X : x value.
C           Q2 : Q2 value.
C           P2 : P2 value; should be = 0. for an on-shell photon.
C           IP2 : scheme used to evaluate off-shell anomalous component.
C               = 0 : recommended default, see = 7.
C               = 1 : dipole dampening by integration; very time-consuming.
C               = 2 : P_0^2 = max( Q_0^2, P^2 )
C               = 3 : P'_0^2 = Q_0^2 + P^2.
C               = 4 : P_{eff} that preserves momentum sum.
C               = 5 : P_{int} that preserves momentum and average
C                     evolution range.
C               = 6 : P_{eff}, matched to P_0 in P2 -> Q2 limit.
C               = 7 : P_{eff}, matched to P_0 in P2 -> Q2 limit.
C...Output: F2GM : F2 value of the photon (including factors of alpha_em).
C           XPFDGM :  x times parton distribution functions of the photon,
C               with elements 0 = g, 1 = d, 2 = u, 3 = s, 4 = c, 5 = b,
C               6 = t (always empty!), - for antiquarks (result is same).
C...The breakdown by component is stored in the commonblock SASCOM,
C               with elements as above.
C           XPVMD : rho, omega, phi VMD part only of output.
C           XPANL : d, u, s anomalous part only of output.
C           XPANH : c, b anomalous part only of output.
C           XPBEH : c, b Bethe-Heitler part only of output.
C           XPDIR : Cgamma (direct contribution) part only of output.
C...The above arrays do not distinguish valence and sea contributions,
C...although this information is available internally. The additional
C...commonblock SASVAL provides the valence part only of the above
C...distributions. Array names VXPVMD, VXPANL and VXPANH correspond
C...to XPVMD, XPANL and XPANH, while XPBEH and XPDIR are valence only
C...and therefore not given doubly. VXPDGM gives the sum of valence
C...parts, and so matches XPDFGM. The difference, i.e. XPVMD-VXPVMD
C...and so on, gives the sea part only.
 
      SUBROUTINE SASGAM2(ISET,X,Q2,P2,IP2,F2GM,XPDFGM)
C...Purpose: to construct the F2 and parton distributions of the photon
C...by summing homogeneous (VMD) and inhomogeneous (anomalous) terms.
C...For F2, c and b are included by the Bethe-Heitler formula;
C...in the 'MSbar' scheme additionally a Cgamma term is added.
      DIMENSION XPDFGM(-6:6)
      COMMON/SASCOM/XPVMD(-6:6),XPANL(-6:6),XPANH(-6:6),XPBEH(-6:6),
     &XPDIR(-6:6)
      COMMON/SASVAL/VXPVMD(-6:6),VXPANL(-6:6),VXPANH(-6:6),VXPDGM(-6:6)
      SAVE /SASCOM/,/SASVAL/
 
C...Temporary array.
      DIMENSION XPGA(-6:6), VXPGA(-6:6)
C...Charm and bottom masses (low to compensate for J/psi etc.).
      DATA PMC/1.3/, PMB/4.6/
C...alpha_em and alpha_em/(2*pi).
      DATA AEM/0.007297/, AEM2PI/0.0011614/
C...Lambda value for 4 flavours.
      DATA ALAM/0.20/
C...Mixture u/(u+d), = 0.5 for incoherent and = 0.8 for coherent sum.
      DATA FRACU/0.8/
C...VMD couplings f_V**2/(4*pi).
      DATA FRHO/2.20/, FOMEGA/23.6/, FPHI/18.4/
C...Masses for rho (=omega) and phi.
      DATA PMRHO/0.770/, PMPHI/1.020/
C...Number of points in integration for IP2=1.
      DATA NSTEP/100/
 
C...Reset output.
      F2GM=0.
      DO 100 KFL=-6,6
      XPDFGM(KFL)=0.
      XPVMD(KFL)=0.
      XPANL(KFL)=0.
      XPANH(KFL)=0.
      XPBEH(KFL)=0.
      XPDIR(KFL)=0.
      VXPVMD(KFL)=0.
      VXPANL(KFL)=0.
      VXPANH(KFL)=0.
      VXPDGM(KFL)=0.
  100 CONTINUE
 
C...Check that input sensible.
      IF(ISET.LE.0.OR.ISET.GE.5) THEN
        WRITE(*,*) ' FATAL ERROR: SaSgam called for unknown set'
        WRITE(*,*) ' ISET = ',ISET
        STOP
      ENDIF
      IF(X.LE.0..OR.X.GT.1.) THEN
        WRITE(*,*) ' FATAL ERROR: SaSgam called for unphysical x'
        WRITE(*,*) ' X = ',X
        STOP
      ENDIF
 
C...Set Q0 cut-off parameter as function of set used.
      IF(ISET.LE.2) THEN
        Q0=0.6
      ELSE
        Q0=2.
      ENDIF
      Q02=Q0**2
 
C...Scale choice for off-shell photon; common factors.
      Q2A=Q2
      FACNOR=1.
      IF(IP2.EQ.1) THEN
        P2MX=P2+Q02
        Q2A=Q2+P2*Q02/MAX(Q02,Q2)
        FACNOR=LOG(Q2/Q02)/NSTEP
      ELSEIF(IP2.EQ.2) THEN
        P2MX=MAX(P2,Q02)
      ELSEIF(IP2.EQ.3) THEN
        P2MX=P2+Q02
        Q2A=Q2+P2*Q02/MAX(Q02,Q2)
      ELSEIF(IP2.EQ.4) THEN
        P2MX=Q2*(Q02+P2)/(Q2+P2)*EXP(P2*(Q2-Q02)/
     &  ((Q2+P2)*(Q02+P2)))
      ELSEIF(IP2.EQ.5) THEN
        P2MXA=Q2*(Q02+P2)/(Q2+P2)*EXP(P2*(Q2-Q02)/
     &  ((Q2+P2)*(Q02+P2)))
        P2MX=Q0*SQRT(P2MXA)
        FACNOR=LOG(Q2/P2MXA)/LOG(Q2/P2MX)
      ELSEIF(IP2.EQ.6) THEN
        P2MX=Q2*(Q02+P2)/(Q2+P2)*EXP(P2*(Q2-Q02)/
     &  ((Q2+P2)*(Q02+P2)))
        P2MX=MAX(0.,1.-P2/Q2)*P2MX+MIN(1.,P2/Q2)*MAX(P2,Q02)
      ELSE
        P2MXA=Q2*(Q02+P2)/(Q2+P2)*EXP(P2*(Q2-Q02)/
     &  ((Q2+P2)*(Q02+P2)))
        P2MX=Q0*SQRT(P2MXA)
        P2MXB=P2MX
        P2MX=MAX(0.,1.-P2/Q2)*P2MX+MIN(1.,P2/Q2)*MAX(P2,Q02)
        P2MXB=MAX(0.,1.-P2/Q2)*P2MXB+MIN(1.,P2/Q2)*P2MXA
        FACNOR=LOG(Q2/P2MXA)/LOG(Q2/P2MXB)
      ENDIF
 
C...Call VMD parametrization for d quark and use to give rho, omega,
C...phi. Note dipole dampening for off-shell photon.
      CALL SASVMD(ISET,1,X,Q2A,P2MX,ALAM,XPGA,VXPGA)
      XFVAL=VXPGA(1)
      XPGA(1)=XPGA(2)
      XPGA(-1)=XPGA(-2)
      FACUD=AEM*(1./FRHO+1./FOMEGA)*(PMRHO**2/(PMRHO**2+P2))**2
      FACS=AEM*(1./FPHI)*(PMPHI**2/(PMPHI**2+P2))**2
      DO 110 KFL=-5,5
      XPVMD(KFL)=(FACUD+FACS)*XPGA(KFL)
  110 CONTINUE
      XPVMD(1)=XPVMD(1)+(1.-FRACU)*FACUD*XFVAL
      XPVMD(2)=XPVMD(2)+FRACU*FACUD*XFVAL
      XPVMD(3)=XPVMD(3)+FACS*XFVAL
      XPVMD(-1)=XPVMD(-1)+(1.-FRACU)*FACUD*XFVAL
      XPVMD(-2)=XPVMD(-2)+FRACU*FACUD*XFVAL
      XPVMD(-3)=XPVMD(-3)+FACS*XFVAL
      VXPVMD(1)=(1.-FRACU)*FACUD*XFVAL
      VXPVMD(2)=FRACU*FACUD*XFVAL
      VXPVMD(3)=FACS*XFVAL
      VXPVMD(-1)=(1.-FRACU)*FACUD*XFVAL
      VXPVMD(-2)=FRACU*FACUD*XFVAL
      VXPVMD(-3)=FACS*XFVAL
 
      IF(IP2.NE.1) THEN
C...Anomalous parametrizations for different strategies
C...for off-shell photons; except full integration.
 
C...Call anomalous parametrization for d + u + s.
        CALL SASANO(-3,X,Q2A,P2MX,ALAM,XPGA,VXPGA)
        DO 120 KFL=-5,5
        XPANL(KFL)=FACNOR*XPGA(KFL)
        VXPANL(KFL)=FACNOR*VXPGA(KFL)
  120   CONTINUE
 
C...Call anomalous parametrization for c and b.
        CALL SASANO(4,X,Q2A,P2MX,ALAM,XPGA,VXPGA)
        DO 130 KFL=-5,5
        XPANH(KFL)=FACNOR*XPGA(KFL)
        VXPANH(KFL)=FACNOR*VXPGA(KFL)
  130   CONTINUE
        CALL SASANO(5,X,Q2A,P2MX,ALAM,XPGA,VXPGA)
        DO 140 KFL=-5,5
        XPANH(KFL)=XPANH(KFL)+FACNOR*XPGA(KFL)
        VXPANH(KFL)=VXPANH(KFL)+FACNOR*VXPGA(KFL)
  140   CONTINUE
 
      ELSE
C...Special option: loop over flavours and integrate over k2.
        DO 170 KF=1,5
        DO 160 ISTEP=1,NSTEP
        Q2STEP=Q02*(Q2/Q02)**((ISTEP-0.5)/NSTEP)
        IF((KF.EQ.4.AND.Q2STEP.LT.PMC**2).OR.
     &  (KF.EQ.5.AND.Q2STEP.LT.PMB**2)) GOTO 160
        CALL SASVMD(0,KF,X,Q2,Q2STEP,ALAM,XPGA,VXPGA)
        FACQ=AEM2PI*(Q2STEP/(Q2STEP+P2))**2*FACNOR
        IF(MOD(KF,2).EQ.0) FACQ=FACQ*(8./9.)
        IF(MOD(KF,2).EQ.1) FACQ=FACQ*(2./9.)
        DO 150 KFL=-5,5
        IF(KF.LE.3) XPANL(KFL)=XPANL(KFL)+FACQ*XPGA(KFL)
        IF(KF.GE.4) XPANH(KFL)=XPANH(KFL)+FACQ*XPGA(KFL)
        IF(KF.LE.3) VXPANL(KFL)=VXPANL(KFL)+FACQ*VXPGA(KFL)
        IF(KF.GE.4) VXPANH(KFL)=VXPANH(KFL)+FACQ*VXPGA(KFL)
  150   CONTINUE
  160   CONTINUE
  170   CONTINUE
      ENDIF
 
C...Call Bethe-Heitler term expression for charm and bottom.
      CALL SASBEH(4,X,Q2,P2,PMC**2,XPBH)
      XPBEH(4)=XPBH
      XPBEH(-4)=XPBH
      CALL SASBEH(5,X,Q2,P2,PMB**2,XPBH)
      XPBEH(5)=XPBH
      XPBEH(-5)=XPBH
 
C...For MSbar subtraction call C^gamma term expression for d, u, s.
      IF(ISET.EQ.2.OR.ISET.EQ.4) THEN
        CALL SASDIR(X,Q2,P2,Q02,XPGA)
        DO 180 KFL=-5,5
        XPDIR(KFL)=XPGA(KFL)
  180   CONTINUE
      ENDIF
 
C...Store result in output array.
      DO 190 KFL=-5,5
      CHSQ=1./9.
      IF(IABS(KFL).EQ.2.OR.IABS(KFL).EQ.4) CHSQ=4./9.
      XPF2=XPVMD(KFL)+XPANL(KFL)+XPBEH(KFL)+XPDIR(KFL)
      IF(KFL.NE.0) F2GM=F2GM+CHSQ*XPF2
      XPDFGM(KFL)=XPVMD(KFL)+XPANL(KFL)+XPANH(KFL)
      VXPDGM(KFL)=VXPVMD(KFL)+VXPANL(KFL)+VXPANH(KFL)
  190 CONTINUE
 
      RETURN
      END
