 
C*********************************************************************
 
C...The following routines are adapted from
C...SaSgam - parton distributions of the photon
C...by Gerhard A. Schuler and Torbjorn Sjostrand
C...For further information see CERN-TH/95-62.
C...The version found here is NOT suitable for standalone usage.
 
      SUBROUTINE PYGGAM(ISET,X,Q2,P2,F2GM,XPDFGM)
C...Purpose: to construct the F2 and parton distributions of the photon
C...by summing homogeneous (VMD) and inhomogeneous (anomalous) terms.
C...For F2, c and b are included by the Bethe-Heitler formula;
C...in the 'MSbar' scheme additionally a Cgamma term is added.
      DIMENSION XPDFGM(-6:6)
      COMMON/PYINT8/XPVMD(-6:6),XPANL(-6:6),XPANH(-6:6),XPBEH(-6:6),
     &XPDIR(-6:6)
      SAVE /PYINT8/
 
C...Temporary array.
      DIMENSION XPGA(-6:6)
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
 
C...Reset output.
      F2GM=0.
      DO 100 KFL=-6,6
      XPDFGM(KFL)=0.
      XPVMD(KFL)=0.
      XPANL(KFL)=0.
      XPANH(KFL)=0.
      XPBEH(KFL)=0.
      XPDIR(KFL)=0.
  100 CONTINUE
 
C...Set k0 cut-off parameter as function of set used.
      IF(ISET.LE.2) THEN
        AK0=0.6
      ELSE
        AK0=2.
      ENDIF
 
C...Call VMD parametrization for d quark and use to give rho+omega+ phi.
C...Note scale choice and dipole dampening for off-shell photon.
      P2MX=MAX(P2,AK0**2)
      CALL PYGVMD(ISET,1,X,Q2,P2MX,ALAM,XPGA)
      XFVAL=XPGA(1)-XPGA(2)
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
 
C...Call anomalous parametrization for d + u + s.
      CALL PYGANO(-3,X,Q2,P2MX,ALAM,XPGA)
      DO 120 KFL=-5,5
      XPANL(KFL)=XPGA(KFL)
  120 CONTINUE
 
C...Call anomalous parametrization for c and b.
      CALL PYGANO(4,X,Q2,P2MX,ALAM,XPGA)
      DO 130 KFL=-5,5
      XPANH(KFL)=XPGA(KFL)
  130 CONTINUE
      CALL PYGANO(5,X,Q2,P2MX,ALAM,XPGA)
      DO 140 KFL=-5,5
      XPANH(KFL)=XPANH(KFL)+XPGA(KFL)
  140 CONTINUE
 
C...Call Bethe-Heitler term expression for charm and bottom.
      CALL PYGBEH(4,X,Q2,P2,PMC**2,XPBH)
      XPBEH(4)=XPBH
      XPBEH(-4)=XPBH
      CALL PYGBEH(5,X,Q2,P2,PMB**2,XPBH)
      XPBEH(5)=XPBH
      XPBEH(-5)=XPBH
 
C...For MSbar subtraction call C^gamma term expression for d, u, s.
      IF(ISET.EQ.2.OR.ISET.EQ.4) THEN
        CALL PYGDIR(X,Q2,P2,AK0,XPGA)
        DO 150 KFL=-5,5
        XPDIR(KFL)=XPGA(KFL)
  150   CONTINUE
      ENDIF
 
C...Store result in output array.
      DO 160 KFL=-5,5
      CHSQ=1./9.
      IF(IABS(KFL).EQ.2.OR.IABS(KFL).EQ.4) CHSQ=4./9.
      XPF2=XPVMD(KFL)+XPANL(KFL)+XPBEH(KFL)+XPDIR(KFL)
      IF(KFL.NE.0) F2GM=F2GM+CHSQ*XPF2
      XPDFGM(KFL)=XPVMD(KFL)+XPANL(KFL)+XPANH(KFL)
  160 CONTINUE
 
      RETURN
      END
