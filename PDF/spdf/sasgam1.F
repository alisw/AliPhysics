#include "pdf/pilot.h"
C...SaSgam - parton distributions of the photon
C...by Gerhard A. Schuler and Torbjorn Sjostrand
C...For further information see preprint CERN-TH/95-62 and LU TP 95-6:
C...Low- and high-mass components of the photon distribution functions
C...Program last changed on 21 March 1995.
 
C...The user should only need to call the SASGAM routine,
C...which in turn calls the auxiliary routines SASVM1, SASAN1,
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
C     CALL SASGAM1(ISET,X,Q2,P2,F2GM,XPDFGM)
C...with the following declaration statement:
C     DIMENSION XPDFGM(-6:6)
C...and, optionally, further information in:
C     COMMON/SASCOM/XPVMD(-6:6),XPANL(-6:6),XPANH(-6:6),XPBEH(-6:6),
C    &XPDIR(-6:6)
C...Input:  ISET = 1 : SaS set 1D ('DIS',   Q0 = 0.6 GeV)
C                = 2 : SaS set 1M ('MSbar', Q0 = 0.6 GeV)
C                = 3 : SaS set 2D ('DIS',   Q0 =  2  GeV)
C                = 4 : SaS set 2M ('MSbar', Q0 =  2  GeV)
C           X : x value.
C           Q2 : Q2 value.
C           P2 : P2 value; should be = 0. for an on-shell photon.
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
 
      SUBROUTINE SASGAM1(ISET,X,Q2,P2,F2GM,XPDFGM)
C...Purpose: to construct the F2 and parton distributions of the photon
C...by summing homogeneous (VMD) and inhomogeneous (anomalous) terms.
C...For F2, c and b are included by the Bethe-Heitler formula;
C...in the 'MSbar' scheme additionally a Cgamma term is added.
      DIMENSION XPDFGM(-6:6)
      COMMON/SASCOM/XPVMD(-6:6),XPANL(-6:6),XPANH(-6:6),XPBEH(-6:6),
     &XPDIR(-6:6)
      SAVE /SASCOM/
 
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
      Q02 = Q0**2

C...Call VMD parametrization for d quark and use to give rho, omega, phi.
C...Note scale choice and dipole dampening for off-shell photon.
      P2MX=MAX(P2,Q02)
      CALL SASVM1(ISET,1,X,Q2,P2MX,ALAM,XPGA)
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
      CALL SASAN1(-3,X,Q2,P2MX,ALAM,XPGA)
      DO 120 KFL=-5,5
      XPANL(KFL)=XPGA(KFL)
  120 CONTINUE

C...Call anomalous parametrization for c and b.
      CALL SASAN1(4,X,Q2,P2MX,ALAM,XPGA)
      DO 130 KFL=-5,5
      XPANH(KFL)=XPGA(KFL)
  130 CONTINUE
      CALL SASAN1(5,X,Q2,P2MX,ALAM,XPGA)
      DO 140 KFL=-5,5
      XPANH(KFL)=XPANH(KFL)+XPGA(KFL)
  140 CONTINUE
 
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
