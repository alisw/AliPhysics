      subroutine SASGevolvep(xin,qin,p2in,ip2in,pdf)
      real*8 xin,qin,q2in,p2in,pdf(-6:6),xval(45),qcdl4,qcdl5     
      include 'parmsetup.inc'
      character*16 name(nmxset)
      integer nmem(nmxset),ndef(nmxset),mmem
      common/NAME/name,nmem,ndef,mmem
      integer nset
      real*8 upv,dnv,usea,dsea,str,chm,bot,top,glu
      
      save 
      call getnset(iset)
      call getnmem(iset,imem)
      
      iimem = imem
      if(iimem.eq.0) iimem=6
      q2in = qin*qin
      call SFSASxx(iimem,xin,Q2in,p2in,ip2,
     +                     upv,dnv,usea,dsea,str,chm,bot,top,glu)
      
      pdf(-6)= top
      pdf(6)= top
      pdf(-5)= bot
      pdf(5 )= bot
      pdf(-4)= chm
      pdf(4 )= chm
      pdf(-3)= str
      pdf(3 )= str
      pdf(-2)= usea
      pdf(2 )= upv
      pdf(-1)= dsea
      pdf(1 )= dnv
      pdf(0 )= glu

      return
c      
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      entry SASGread(nset)
c      print *,'calling SASGread'
      read(1,*)nmem(nset),ndef(nset)
      return
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      entry SASGalfa(alfas,qalfa)
        call getnset(iset)
	call getnmem(iset,imem)
	call GetOrderAsM(iset,iord)
        call Getlam4M(iset,imem,qcdl4)
        call Getlam5M(iset,imem,qcdl5)
        call aspdflib(alfas,Qalfa,iord,qcdl5)

      return
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      entry SASGinit(Eorder,Q2fit)
      return
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      entry SASGpdf(mem)
c      print *,'calling SASGpdf',mem
c      imem = mem
      call getnset(iset)
      call setnmem(iset,mem)
      return
c
 1000 format(5e13.5)
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C-----------------------------------------------------------------------
      SUBROUTINE SFSASxx(iset,DX,DQ2,DP2,ip2,
     +                     DUPV,DDNV,DSEA,DSEAD,DSTR,DCHM,DBOT,DTOP,DGL)
C
C   ********************************************************************
C   *                                                                  *
C   *        Interface to SASset of structure functions                *
C   *                                                                  *
C   *        Author:    H. Plothow-Besch (CERN-PPE)                    *
C   *                                                                  *
C   ********************************************************************
C
C :::::::::::: Structure functions from the SAS group version 2
C :::::::::::: Lambda = 0.200 GeV, Q**2 = 0.36 GeV**2 (DIS)
C
      double precision
     +          DX,DQ2,DP2,
     +          DUPV,DDNV,DSEA,DSEAD,DSTR,DCHM,DBOT,DTOP,DGL
      DIMENSION XPDFGM(-6:6)
      REAL X, Q, Q2, P2, F2GAM, XPDFGM
c      PARAMETER (ISET=1)
C
      X  = DX
      Q  = SQRT(DQ2)
      Q2 = DQ2
      P2 = DP2
C
C     generate the individual structure fcn calls
C
      if(iset.le.4) then
         iiset=iset
         CALL LHASASGAM1(iISET,X,Q2,P2,F2GAM,XPDFGM)
      else
         iiset = iset-4
         CALL LHASASGAM2(iISET,X,Q2,P2,ip2,F2GAM,XPDFGM)
      endif
      
      UPV = XPDFGM(2)
      DUPV = UPV
      DNV = XPDFGM(1)
      DDNV = DNV
      SEAU = XPDFGM(-2)
      DSEA = SEAU
      SEAD = XPDFGM(-1)
      DSEAD = SEAD
      STR = XPDFGM(-3)
      DSTR = STR
      CHM = XPDFGM(-4)
      DCHM = CHM
      BOT = XPDFGM(-5)
      DBOT = BOT
      TOP = 0.
C      IF (DSCAL.GT.TMAS) TOP = XPDFGM(6)
      DTOP = TOP
      GL = XPDFGM(0)
      DGL = GL
C
      RETURN
      END
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c ------------- SASGAM1 ------------------------
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
 
      SUBROUTINE LHASASGAM1(ISET,X,Q2,P2,F2GM,XPDFGM)
C...Purpose: to construct the F2 and parton distributions of the photon
C...by summing homogeneous (VMD) and inhomogeneous (anomalous) terms.
C...For F2, c and b are included by the Bethe-Heitler formula;
C...in the 'MSbar' scheme additionally a Cgamma term is added.
      DIMENSION XPDFGM(-6:6)
      COMMON/LHASASCOM/XPVMD(-6:6),XPANL(-6:6),XPANH(-6:6),XPBEH(-6:6),
     &XPDIR(-6:6)
      SAVE /LHASASCOM/
 
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
      CALL LHASASVM1(ISET,1,X,Q2,P2MX,ALAM,XPGA)
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
      CALL LHASASAN1(-3,X,Q2,P2MX,ALAM,XPGA)
      DO 120 KFL=-5,5
      XPANL(KFL)=XPGA(KFL)
  120 CONTINUE

C...Call anomalous parametrization for c and b.
      CALL LHASASAN1(4,X,Q2,P2MX,ALAM,XPGA)
      DO 130 KFL=-5,5
      XPANH(KFL)=XPGA(KFL)
  130 CONTINUE
      CALL LHASASAN1(5,X,Q2,P2MX,ALAM,XPGA)
      DO 140 KFL=-5,5
      XPANH(KFL)=XPANH(KFL)+XPGA(KFL)
  140 CONTINUE
 
C...Call Bethe-Heitler term expression for charm and bottom.
      CALL LHASASBEH(4,X,Q2,P2,PMC**2,XPBH)
      XPBEH(4)=XPBH
      XPBEH(-4)=XPBH
      CALL LHASASBEH(5,X,Q2,P2,PMB**2,XPBH)
      XPBEH(5)=XPBH
      XPBEH(-5)=XPBH

C...For MSbar subtraction call C^gamma term expression for d, u, s.
      IF(ISET.EQ.2.OR.ISET.EQ.4) THEN
        CALL LHASASDIR(X,Q2,P2,Q02,XPGA)
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
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c ------------------ SASGAM2 ---------------------------
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
 
      SUBROUTINE LHASASGAM2(ISET,X,Q2,P2,IP2,F2GM,XPDFGM)
C...Purpose: to construct the F2 and parton distributions of the photon
C...by summing homogeneous (VMD) and inhomogeneous (anomalous) terms.
C...For F2, c and b are included by the Bethe-Heitler formula;
C...in the 'MSbar' scheme additionally a Cgamma term is added.
      DIMENSION XPDFGM(-6:6)
      COMMON/LHASASCOM/XPVMD(-6:6),XPANL(-6:6),XPANH(-6:6),XPBEH(-6:6),
     &XPDIR(-6:6)
      COMMON/LHASASVAL/VXPVMD(-6:6),VXPANL(-6:6),VXPANH(-6:6),
     &VXPDGM(-6:6)
      SAVE /LHASASCOM/,/LHASASVAL/
 
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
      CALL LHASASVMD(ISET,1,X,Q2A,P2MX,ALAM,XPGA,VXPGA)
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
        CALL LHASASANO(-3,X,Q2A,P2MX,ALAM,XPGA,VXPGA)
        DO 120 KFL=-5,5
        XPANL(KFL)=FACNOR*XPGA(KFL)
        VXPANL(KFL)=FACNOR*VXPGA(KFL)
  120   CONTINUE
 
C...Call anomalous parametrization for c and b.
        CALL LHASASANO(4,X,Q2A,P2MX,ALAM,XPGA,VXPGA)
        DO 130 KFL=-5,5
        XPANH(KFL)=FACNOR*XPGA(KFL)
        VXPANH(KFL)=FACNOR*VXPGA(KFL)
  130   CONTINUE
        CALL LHASASANO(5,X,Q2A,P2MX,ALAM,XPGA,VXPGA)
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
        CALL LHASASVMD(0,KF,X,Q2,Q2STEP,ALAM,XPGA,VXPGA)
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
      CALL LHASASBEH(4,X,Q2,P2,PMC**2,XPBH)
      XPBEH(4)=XPBH
      XPBEH(-4)=XPBH
      CALL LHASASBEH(5,X,Q2,P2,PMB**2,XPBH)
      XPBEH(5)=XPBH
      XPBEH(-5)=XPBH
 
C...For MSbar subtraction call C^gamma term expression for d, u, s.
      IF(ISET.EQ.2.OR.ISET.EQ.4) THEN
        CALL LHASASDIR(X,Q2,P2,Q02,XPGA)
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


ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      SUBROUTINE LHASASVMD(ISET,KF,X,Q2,P2,ALAM,XPGA,VXPGA)
C...Purpose: to evaluate the VMD parton distributions of a photon,
C...evolved homogeneously from an initial scale P2 to Q2.
C...Does not include dipole suppression factor.
C...ISET is parton distribution set, see above; 
C...additionally ISET=0 is used for the evolution of an anomalous photon 
C...which branched at a scale P2 and then evolved homogeneously to Q2.
C...ALAM is the 4-flavour Lambda, which is automatically converted
C...to 3- and 5-flavour equivalents as needed.
      DIMENSION XPGA(-6:6), VXPGA(-6:6)
      DATA PMC/1.3/, PMB/4.6/, AEM/0.007297/, AEM2PI/0.0011614/

C...Reset output.
      DO 100 KFL=-6,6
      XPGA(KFL)=0.
      VXPGA(KFL)=0.
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
      VXPGA(KFA)=XVAL
      VXPGA(-KFA)=XVAL

      RETURN
      END 
 
C********************************************************************* 

      SUBROUTINE LHASASANO(KF,X,Q2,P2,ALAM,XPGA,VXPGA)
C...Purpose: to evaluate the parton distributions of the anomalous 
C...photon, inhomogeneously evolved from a scale P2 (where it vanishes) 
C...to Q2.
C...KF=0 gives the sum over (up to) 5 flavours,
C...KF<0 limits to flavours up to abs(KF),
C...KF>0 is for flavour KF only. 
C...ALAM is the 4-flavour Lambda, which is automatically converted
C...to 3- and 5-flavour equivalents as needed.
      DIMENSION XPGA(-6:6), VXPGA(-6:6), ALAMSQ(3:5)
      DATA PMC/1.3/, PMB/4.6/, AEM/0.007297/, AEM2PI/0.0011614/

C...Reset output.
      DO 100 KFL=-6,6
      XPGA(KFL)=0.
      VXPGA(KFL)=0.
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
      VXPGA(KFL)=VXPGA(KFL)+FAC*XVAL
  110 CONTINUE
      DO 120 KFL=1,5
      XPGA(-KFL)=XPGA(KFL)
      VXPGA(-KFL)=VXPGA(KFL)
  120 CONTINUE

      RETURN
      END 
 
C********************************************************************* 

      SUBROUTINE LHASASBEH(KF,X,Q2,P2,PM2,XPBH)
C...Purpose: to evaluate the Bethe-Heitler cross section for
C...heavy flavour production.
      DATA AEM2PI/0.0011614/

C...Reset output.
      XPBH=0.
      SIGBH=0.

C...Check kinematics limits.
      IF(X.GE.Q2/(4.*PM2+Q2+P2)) RETURN           
      W2=Q2*(1.-X)/X-P2
      BETA2=1.-4.*PM2/W2
      IF(BETA2.LT.1E-10) RETURN
      BETA=SQRT(BETA2)
      RMQ=4.*PM2/Q2
 
C...Simple case: P2 = 0.
      IF(P2.LT.1E-4) THEN
        IF(BETA.LT.0.99) THEN
          XBL=LOG((1.+BETA)/(1.-BETA))
        ELSE
          XBL=LOG((1.+BETA)**2*W2/(4.*PM2))
        ENDIF 
        SIGBH=BETA*(8.*X*(1.-X)-1.-RMQ*X*(1.-X))+
     &  XBL*(X**2+(1.-X)**2+RMQ*X*(1.-3.*X)-0.5*RMQ**2*X**2)
    
C...Complicated case: P2 > 0, based on approximation of
C...C.T. Hill and G.G. Ross, Nucl. Phys. B148 (1979) 373
      ELSE
        RPQ=1.-4.*X**2*P2/Q2
        IF(RPQ.GT.1E-10) THEN
          RPBE=SQRT(RPQ*BETA2)
          IF(RPBE.LT.0.99) THEN
            XBL=LOG((1.+RPBE)/(1.-RPBE))
            XBI=2.*RPBE/(1.-RPBE**2)
          ELSE
            RPBESN=4.*PM2/W2+(4.*X**2*P2/Q2)*BETA2
            XBL=LOG((1.+RPBE)**2/RPBESN)
            XBI=2.*RPBE/RPBESN
          ENDIF
          SIGBH=BETA*(6.*X*(1.-X)-1.)+
     &    XBL*(X**2+(1.-X)**2+RMQ*X*(1.-3.*X)-0.5*RMQ**2*X**2)+
     &    XBI*(2.*X/Q2)*(PM2*X*(2.-RMQ)-P2*X)
        ENDIF                
      ENDIF

C...Multiply by charge-squared etc. to get parton distribution.
      CHSQ=1./9.
      IF(IABS(KF).EQ.2.OR.IABS(KF).EQ.4) CHSQ=4./9.
      XPBH=3.*CHSQ*AEM2PI*X*SIGBH       

      RETURN
      END
 
C********************************************************************* 

       SUBROUTINE LHASASDIR(X,Q2,P2,Q02,XPGA)
C...Purpose: to evaluate the direct contribution, i.e. the C^gamma term,
C...as needed in MSbar parametrizations.
      DIMENSION XPGA(-6:6)
      DATA PMC/1.3/, PMB/4.6/, AEM2PI/0.0011614/

C...Reset output.
      DO 100 KFL=-6,6
      XPGA(KFL)=0.
  100 CONTINUE

C...Evaluate common x-dependent expression.
      XTMP = (X**2+(1.-X)**2) * (-LOG(X)) - 1.
      CGAM = 3.*AEM2PI*X * (XTMP*(1.+P2/(P2+Q02)) + 6.*X*(1.-X))

C...d, u, s part by simple charge factor.
      XPGA(1)=(1./9.)*CGAM
      XPGA(2)=(4./9.)*CGAM
      XPGA(3)=(1./9.)*CGAM      

C...Also fill for antiquarks.     
      DO 110 KF=1,5
      XPGA(-KF)=XPGA(KF)
  110 CONTINUE

      RETURN
      END
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
       SUBROUTINE LHASASVM1(ISET,KF,X,Q2,P2,ALAM,XPGA)
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
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      SUBROUTINE LHASASAN1(KF,X,Q2,P2,ALAM,XPGA)
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
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
