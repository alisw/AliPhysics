#include "isajet/pilot.h"
C  *********************************************************************
C  ***                                                               ***
C  ***               coded by H. Murayama & I. Watanabe              ***
C  *** For the formalism and notations, see the following reference: ***
C  ***           H. Murayama, I. Watanabe and K. Hagiwara            ***
C  ***           "HELAS: HELicity Amplitude Subroutines              ***
C  ***               for Feynman diagram evaluation"                 ***
C  ***               KEK Report 91-11, December 1991                 ***
C  ***                                                               ***
C  *********************************************************************
C
C  Converted to double precision by W. Long and T. Seltzer for MadGraph.
C
C  Minor changes for portability by FEP, July 1999. The code is not ANSI
C  standard, but that cannot be helped if MadGraph compatibility is to 
C  be maintained.
C
C ======================================================================
C
      SUBROUTINE BOOSTX(P,Q , PBOOST)
C
C this subroutine performs the lorentz boost of a four-momentum.  the   
C momentum p is assumed to be given in the rest frame of q.  pboost is  
C the momentum p boosted to the frame in which q is given.  q must be a 
C timelike momentum.                                                    
C                                                                       
C input:                                                                
C       real    p(0:3)         : four-momentum p in the q rest  frame   
C       real    q(0:3)         : four-momentum q in the boosted frame   
C                                                                       
C output:                                                               
C       real    pboost(0:3)    : four-momentum p in the boosted frame   
C
#if defined(CERNLIB_IMPNONE)
      IMPLICIT NONE
#endif
      REAL*8    P(0:3),Q(0:3),PBOOST(0:3),PQ,QQ,M,LF
      REAL*8 RXZERO
      PARAMETER( RXZERO=0.0D0 )
C
      QQ=Q(1)**2+Q(2)**2+Q(3)**2
C
      IF ( QQ .NE. RXZERO ) THEN
         PQ=P(1)*Q(1)+P(2)*Q(2)+P(3)*Q(3)
         M=SQRT(Q(0)**2-QQ)
         LF=((Q(0)-M)*PQ/QQ+P(0))/M
         PBOOST(0) = (P(0)*Q(0)+PQ)/M
         PBOOST(1) =  P(1)+Q(1)*LF
         PBOOST(2) =  P(2)+Q(2)*LF
         PBOOST(3) =  P(3)+Q(3)*LF
      ELSE
         PBOOST(0)=P(0)
         PBOOST(1)=P(1)
         PBOOST(2)=P(2)
         PBOOST(3)=P(3)
      ENDIF
C
      RETURN
      END
C
C **********************************************************************
C
      SUBROUTINE COUP1X(SW2 , GW,GWWA,GWWZ)
C
C this subroutine sets up the coupling constants of the gauge bosons in 
C the standard model.                                                   
C                                                                       
C input:                                                                
C       real    sw2            : square of sine of the weak angle       
C                                                                       
C output:                                                               
C       real    gw             : weak coupling constant                 
C       real    gwwa           : dimensionless coupling of w-,w+,a      
C       real    gwwz           : dimensionless coupling of w-,w+,z      
C
#if defined(CERNLIB_IMPNONE)
      IMPLICIT NONE
#endif
      REAL*8    SW2,GW,GWWA,GWWZ,ALPHA,FOURPI,EE,SW,CW
      REAL*8 RXONE, RXFOUR, RXOTE, RXPI, RIALPH
      PARAMETER( RXONE=1.0D0, RXFOUR=4.0D0, RXOTE=128.0D0 )
      PARAMETER( RXPI=3.14159265358979323846D0, RIALPH=137.0359895D0 )
C
      ALPHA = RXONE / RXOTE
C      alpha = r_one / r_ialph
      FOURPI = RXFOUR * RXPI
      EE=SQRT( ALPHA * FOURPI )
      SW=SQRT( SW2 )
      CW=SQRT( RXONE - SW2 )
C
      GW    =  EE/SW
      GWWA  =  EE
      GWWZ  =  EE*CW/SW
C
      RETURN
      END
C
C ----------------------------------------------------------------------
C
      SUBROUTINE COUP2X(SW2 , GAL,GAU,GAD,GWF,GZN,GZL,GZU,GZD,G1)
C
C this subroutine sets up the coupling constants for the fermion-       
C fermion-vector vertices in the standard model.  the array of the      
C couplings specifies the chirality of the flowing-in fermion.  g??(1)  
C denotes a left-handed coupling, and g??(2) a right-handed coupling.   
C                                                                       
C input:                                                                
C       real    sw2            : square of sine of the weak angle       
C                                                                       
C output:                                                               
C       real    gal(2)         : coupling with a of charged leptons     
C       real    gau(2)         : coupling with a of up-type quarks      
C       real    gad(2)         : coupling with a of down-type quarks    
C       real    gwf(2)         : coupling with w-,w+ of fermions        
C       real    gzn(2)         : coupling with z of neutrinos           
C       real    gzl(2)         : coupling with z of charged leptons     
C       real    gzu(2)         : coupling with z of up-type quarks      
C       real    gzd(2)         : coupling with z of down-type quarks    
C       real    g1(2)          : unit coupling of fermions              
C
#if defined(CERNLIB_IMPNONE)
      IMPLICIT NONE
#endif
      REAL*8 GAL(2),GAU(2),GAD(2),GWF(2),GZN(2),GZL(2),GZU(2),GZD(2),
     &     G1(2),SW2,ALPHA,FOURPI,EE,SW,CW,EZ,EY
C
      REAL*8 RXZERO, RXHALF, RXONE, RXTWO, RTHREE, RXFOUR, RXOTE
      REAL*8 RXPI, RIALPH
      PARAMETER( RXZERO=0.0D0, RXHALF=0.5D0, RXONE=1.0D0, RXTWO=2.0D0,
     $     RTHREE=3.0D0 )
      PARAMETER( RXFOUR=4.0D0, RXOTE=128.0D0 )
      PARAMETER( RXPI=3.14159265358979323846D0, RIALPH=137.0359895D0 )
C
      ALPHA = RXONE / RXOTE
C      alpha = r_one / r_ialph
      FOURPI = RXFOUR * RXPI
      EE=SQRT( ALPHA * FOURPI )
      SW=SQRT( SW2 )
      CW=SQRT( RXONE - SW2 )
      EZ=EE/(SW*CW)
      EY=EE*(SW/CW)
C
      GAL(1) =  EE
      GAL(2) =  EE
      GAU(1) = -EE*RXTWO/RTHREE
      GAU(2) = -EE*RXTWO/RTHREE
      GAD(1) =  EE   /RTHREE
      GAD(2) =  EE   /RTHREE
      GWF(1) = -EE/SQRT(RXTWO*SW2)
      GWF(2) =  RXZERO
      GZN(1) = -EZ*  RXHALF
      GZN(2) =  RXZERO
      GZL(1) = -EZ*(-RXHALF+SW2)
      GZL(2) = -EY
      GZU(1) = -EZ*( RXHALF-SW2*RXTWO/RTHREE)
      GZU(2) =  EY*          RXTWO/RTHREE
      GZD(1) = -EZ*(-RXHALF+SW2   /RTHREE)
      GZD(2) = -EY             /RTHREE
      G1(1)  =  RXONE
      G1(2)  =  RXONE
C
      RETURN
      END
C
C ----------------------------------------------------------------------
C
      SUBROUTINE COUP3X(SW2,ZMASS,HMASS , 
     &                  GWWH,GZZH,GHHH,GWWHH,GZZHH,GHHHH)
C
C this subroutine sets up the coupling constants of the gauge bosons and
C higgs boson in the standard model.                                    
C                                                                       
C input:                                                                
C       real    sw2            : square of sine of the weak angle       
C       real    zmass          : mass of z                              
C       real    hmass          : mass of higgs                          
C                                                                       
C output:                                                               
C       real    gwwh           : dimensionful  coupling of w-,w+,h      
C       real    gzzh           : dimensionful  coupling of z, z, h      
C       real    ghhh           : dimensionful  coupling of h, h, h      
C       real    gwwhh          : dimensionful  coupling of w-,w+,h, h   
C       real    gzzhh          : dimensionful  coupling of z, z, h, h   
C       real    ghhhh          : dimensionless coupling of h, h, h, h   
C
#if defined(CERNLIB_IMPNONE)
      IMPLICIT NONE
#endif
      REAL*8    SW2,ZMASS,HMASS,GWWH,GZZH,GHHH,GWWHH,GZZHH,GHHHH,
     &        ALPHA,FOURPI,EE2,SC2,V
C
      REAL*8 RXHALF, RXONE, RXTWO, RTHREE, RXFOUR, RXOTE
      REAL*8 RXPI, RIALPH
      PARAMETER( RXHALF=0.5D0, RXONE=1.0D0, RXTWO=2.0D0, RTHREE=3.0D0 )
      PARAMETER( RXFOUR=4.0D0, RXOTE=128.0D0 )
      PARAMETER( RXPI=3.14159265358979323846D0, RIALPH=137.0359895D0 )
C
      ALPHA = RXONE / RXOTE
C      alpha = r_one / r_ialph
      FOURPI = RXFOUR * RXPI
      EE2=ALPHA*FOURPI
      SC2=SW2*( RXONE - SW2 )
      V = RXTWO * ZMASS*SQRT(SC2)/SQRT(EE2)
C
      GWWH  =   EE2/SW2*RXHALF*V
      GZZH  =   EE2/SC2*RXHALF*V
      GHHH  =  -HMASS**2/V*RTHREE
      GWWHH =   EE2/SW2*RXHALF
      GZZHH =   EE2/SC2*RXHALF
      GHHHH = -(HMASS/V)**2*RTHREE
C
      RETURN
      END
C
C ----------------------------------------------------------------------
C
      SUBROUTINE COUP4X(SW2,ZMASS,FMASS , GCHF)
C
C This subroutine sets up the coupling constant for the fermion-fermion-
C Higgs vertex in the STANDARD MODEL.  The coupling is COMPLEX and the  
C array of the coupling specifies the chirality of the flowing-IN       
C fermion.  GCHF(1) denotes a left-handed coupling, and GCHF(2) a right-
C handed coupling.                                                      
C                                                                       
C INPUT:                                                                
C       real    SW2            : square of sine of the weak angle       
C       real    ZMASS          : Z       mass                           
C       real    FMASS          : fermion mass                           
C                                                                       
C OUTPUT:                                                               
C       complex GCHF(2)        : coupling of fermion and Higgs          
C
#if defined(CERNLIB_IMPNONE)
      IMPLICIT NONE
#endif
      COMPLEX*16 GCHF(2)
      REAL*8    SW2,ZMASS,FMASS,ALPHA,FOURPI,EZ,G
C
      ALPHA=1.D0/128.D0
C      ALPHA=1./REAL(137.0359895)
      FOURPI=4.D0*3.14159265358979323846D0
      EZ=SQRT(ALPHA*FOURPI)/SQRT(SW2*(1.D0-SW2))
      G=EZ*FMASS*0.5D0/ZMASS
C
      GCHF(1) = DCMPLX( -G )
      GCHF(2) = DCMPLX( -G )
C
      RETURN
      END
C
C ======================================================================
C
      SUBROUTINE EAIXXX(EB,EA,SHLF,CHLF,PHI,NHE,NHA , EAI)
C
C This subroutine computes an off-shell electron wavefunction after     
C emitting a photon from the electron beam, with a special care for the 
C small angle region.  The momenta are measured in the laboratory frame,
C where the e- beam is along the positive z axis.                       
C                                                                       
C INPUT:                                                                
C       real    EB             : energy (GeV)    of beam  e-            
C       real    EA             : energy (GeV)    of final photon        
C       real    SHLF           : sin(theta/2)    of final photon        
C       real    CHLF           : cos(theta/2)    of final photon        
C       real    PHI            : azimuthal angle of final photon        
C       integer NHE  = -1 or 1 : helicity        of beam  e-            
C       integer NHA  = -1 or 1 : helicity        of final photon        
C                                                                       
C OUTPUT:                                                               
C       complex EAI(6)         : off-shell electron             |e',A,e>
C
#if defined(CERNLIB_IMPNONE)
      IMPLICIT NONE
#endif
      COMPLEX*16 EAI(6),PHS
      REAL*8  EB,EA,SHLF,CHLF,PHI,ME,ALPHA,GAL,RNHE,X,C,S,D,COEFF,
     &        XNNP,XNNM,SNP,CSP
      INTEGER NHE,NHA,NN
C
      ME   = 0.51099906D-3
      ALPHA=1./128.
      GAL  =SQRT(ALPHA*4.*3.14159265D0)
C
      NN=NHA*NHE
      RNHE=NHE
      X=EA/EB
      C=(CHLF+SHLF)*(CHLF-SHLF)
      S=2.*CHLF*SHLF
      D=-1./(EA*EB*(4.*SHLF**2+(ME/EB)**2*C))
      COEFF=-NN*GAL*SQRT(EB)*D
      XNNP=X*(1+NN)
      XNNM=X*(1-NN)
      SNP=SIN(PHI)
      CSP=COS(PHI)
      PHS=DCMPLX( CSP , RNHE*SNP )
C
      EAI((5-3*NHE)/2) = -RNHE*COEFF*ME*S*(1.+XNNP*.5)
      EAI((5-NHE)/2)   =  XNNP*COEFF*ME*CHLF**2*PHS
      EAI((5+NHE)/2)   =  RNHE*COEFF*EB*S*(-2.+XNNM)
      EAI((5+3*NHE)/2) =  XNNM*COEFF*EB*SHLF**2*PHS*2.
C
      EAI(5) =  EB*DCMPLX( 1.-X , 1.-X*C )
      EAI(6) = -EB*X*S*DCMPLX( CSP , SNP )
C
      RETURN
      END
C
C ----------------------------------------------------------------------
C
      SUBROUTINE EAOXXX(EB,EA,SHLF,CHLF,PHI,NHE,NHA , EAO)
C
C This subroutine computes an off-shell positron wavefunction after     
C emitting a photon from the positron beam, with a special care for the 
C small angle region.  The momenta are measured in the laboratory frame,
C where the e+ beam is along the negative z axis.                       
C                                                                       
C INPUT:                                                                
C       real    EB             : energy (GeV)    of beam  e+            
C       real    EA             : energy (GeV)    of final photon        
C       real    SHLF           : sin(theta/2)    of final photon        
C       real    CHLF           : cos(theta/2)    of final photon        
C       real    PHI            : azimuthal angle of final photon        
C       integer NHE  = -1 or 1 : helicity        of beam  e+            
C       integer NHA  = -1 or 1 : helicity        of final photon        
C                                                                       
C OUTPUT:                                                               
C       complex EAO(6)         : off-shell positron             <e,A,e'|
C
#if defined(CERNLIB_IMPNONE)
      IMPLICIT NONE
#endif
      COMPLEX*16 EAO(6),PHS
      REAL*8  EB,EA,SHLF,CHLF,PHI,ME,ALPHA,GAL,RNHE,X,C,S,D,COEFF,
     &        XNNP,XNNM,SNP,CSP
      INTEGER NHE,NHA,NN
C
      ME   = 0.51099906D-3
      ALPHA=1./128.
      GAL  =SQRT(ALPHA*4.*3.14159265D0)
C
      NN=NHA*NHE
      RNHE=NHE
      X=EA/EB
      C=(CHLF+SHLF)*(CHLF-SHLF)
      S=2.*CHLF*SHLF
      D=-1./(EA*EB*(4.*CHLF**2-(ME/EB)**2*C))
      COEFF=NN*GAL*SQRT(EB)*D
      XNNP=X*(1+NN)
      XNNM=X*(1-NN)
      SNP=SIN(PHI)
      CSP=COS(PHI)
      PHS=DCMPLX( CSP ,-RNHE*SNP )
C
      EAO((5-3*NHE)/2) =              COEFF*ME*S*(1.+XNNP*.5)
      EAO((5-NHE)/2)   = RNHE*XNNP    *COEFF*ME*SHLF**2*PHS
      EAO((5+NHE)/2)   =              COEFF*EB*S*(-2.+XNNM)
      EAO((5+3*NHE)/2) = REAL(NHA-NHE)*COEFF*EB*X*CHLF**2*PHS*2.
C
      EAO(5) = EB*DCMPLX( X-1. , X*C+1. )
      EAO(6) = EB*X*S*DCMPLX( CSP , SNP )
C
      RETURN
      END
C
C ----------------------------------------------------------------------
C
      SUBROUTINE FSIXXX(FI,SC,GC,FMASS,FWIDTH , FSI)
C
C this subroutine computes an off-shell fermion wavefunction from a     
C flowing-in external fermion and a vector boson.                       
C                                                                       
C input:                                                                
C       complex*16 fi(6)          : flow-in  fermion                   |fi>
C       complex*16 sc(3)          : input    scalar                      s 
C       complex*16 gc(2)          : coupling constants                 gchf
C       real*8    fmass          : mass  of output fermion f'             
C       real*8    fwidth         : width of output fermion f'             
C                                                                       
C output:                                                               
C       complex fsi(6)         : off-shell fermion             |f',s,fi>
C
#if defined(CERNLIB_IMPNONE)
      IMPLICIT NONE
#endif
      COMPLEX*16 FI(6),SC(3),FSI(6),GC(2),SL1,SL2,SR1,SR2,DS
      REAL*8     PF(0:3),FMASS,FWIDTH,PF2,P0P3,P0M3
C
      FSI(5) = FI(5)-SC(2)
      FSI(6) = FI(6)-SC(3)
C
      PF(0)=DBLE( FSI(5))
      PF(1)=DBLE( FSI(6))
      PF(2)=DIMAG(FSI(6))
      PF(3)=DIMAG(FSI(5))
      PF2=PF(0)**2-(PF(1)**2+PF(2)**2+PF(3)**2)
C
      DS=-SC(1)/DCMPLX(PF2-FMASS**2,MAX(DSIGN(FMASS*FWIDTH ,PF2),0D0))
      P0P3=PF(0)+PF(3)
      P0M3=PF(0)-PF(3)
      SL1=GC(1)*(P0P3*FI(1)+DCONJG(FSI(6))*FI(2))
      SL2=GC(1)*(P0M3*FI(2)      +FSI(6) *FI(1))
      SR1=GC(2)*(P0M3*FI(3)-DCONJG(FSI(6))*FI(4))
      SR2=GC(2)*(P0P3*FI(4)      -FSI(6) *FI(3))
C
      FSI(1) = ( GC(1)*FMASS*FI(1) + SR1 )*DS
      FSI(2) = ( GC(1)*FMASS*FI(2) + SR2 )*DS
      FSI(3) = ( GC(2)*FMASS*FI(3) + SL1 )*DS
      FSI(4) = ( GC(2)*FMASS*FI(4) + SL2 )*DS
C
      RETURN          
      END
C
C ----------------------------------------------------------------------
C
      SUBROUTINE FSOXXX(FO,SC,GC,FMASS,FWIDTH , FSO)
C
C this subroutine computes an off-shell fermion wavefunction from a     
C flowing-out external fermion and a vector boson.                      
C                                                                       
C input:                                                                
C       complex*16 fo(6)          : flow-out fermion                   <fo|
C       complex*16 sc(6)          : input    scalar                      s 
C       complex*16 gc(2)          : coupling constants                 gchf
C       real*8     fmass          : mass  of output fermion f'             
C       real*8     fwidth         : width of output fermion f'             
C                                                                       
C output:                                                               
C       complex fso(6)         : off-shell fermion             <fo,s,f'|
C
#if defined(CERNLIB_IMPNONE)
      IMPLICIT NONE
#endif
      COMPLEX*16 FO(6),SC(6),FSO(6),GC(2),SL1,SL2,SR1,SR2,DS
      REAL*8     PF(0:3),FMASS,FWIDTH,PF2,P0P3,P0M3
C
      FSO(5) = FO(5)+SC(2)
      FSO(6) = FO(6)+SC(3)
C
      PF(0)=DBLE( FSO(5))
      PF(1)=DBLE( FSO(6))
      PF(2)=DIMAG(FSO(6))
      PF(3)=DIMAG(FSO(5))
      PF2=PF(0)**2-(PF(1)**2+PF(2)**2+PF(3)**2)
C
      DS=-SC(1)/DCMPLX(PF2-FMASS**2,MAX(DSIGN(FMASS*FWIDTH ,PF2),0D0))
      P0P3=PF(0)+PF(3)
      P0M3=PF(0)-PF(3)
      SL1=GC(2)*(P0P3*FO(3)      +FSO(6) *FO(4))
      SL2=GC(2)*(P0M3*FO(4)+DCONJG(FSO(6))*FO(3))
      SR1=GC(1)*(P0M3*FO(1)      -FSO(6) *FO(2))
      SR2=GC(1)*(P0P3*FO(2)-DCONJG(FSO(6))*FO(1))
C
      FSO(1) = ( GC(1)*FMASS*FO(1) + SL1 )*DS
      FSO(2) = ( GC(1)*FMASS*FO(2) + SL2 )*DS
      FSO(3) = ( GC(2)*FMASS*FO(3) + SR1 )*DS
      FSO(4) = ( GC(2)*FMASS*FO(4) + SR2 )*DS
C
      RETURN          
      END
C
C ----------------------------------------------------------------------
C
      SUBROUTINE FVIXXX(FI,VC,G,FMASS,FWIDTH , FVI)
C
C this subroutine computes an off-shell fermion wavefunction from a     
C flowing-in external fermion and a vector boson.                       
C                                                                       
C input:                                                                
C       complex fi(6)          : flow-in  fermion                   |fi>
C       complex vc(6)          : input    vector                      v 
C       real    g(2)           : coupling constants                  gvf
C       real    fmass          : mass  of output fermion f'             
C       real    fwidth         : width of output fermion f'             
C                                                                       
C output:                                                               
C       complex fvi(6)         : off-shell fermion             |f',v,fi>
C
#if defined(CERNLIB_IMPNONE)
      IMPLICIT NONE
#endif
      COMPLEX*16 FI(6),VC(6),FVI(6),SL1,SL2,SR1,SR2,D
      REAL*8    G(2),PF(0:3),FMASS,FWIDTH,PF2
      REAL*8 RXZERO, RXONE
      PARAMETER( RXZERO=0.0D0, RXONE=1.0D0 )
      COMPLEX*16 CXIMAG
C
      LOGICAL FIRST
      SAVE CXIMAG,FIRST
      DATA FIRST/.TRUE./
C
C          Fix compilation with g77
      IF(FIRST) THEN
        FIRST=.FALSE.
        CXIMAG=DCMPLX( RXZERO, RXONE )
      ENDIF
C
      FVI(5) = FI(5)-VC(5)
      FVI(6) = FI(6)-VC(6)
C
      PF(0)=DBLE( FVI(5))
      PF(1)=DBLE( FVI(6))
      PF(2)=DIMAG(FVI(6))
      PF(3)=DIMAG(FVI(5))
      PF2=PF(0)**2-(PF(1)**2+PF(2)**2+PF(3)**2)
C
      D=-RXONE/DCMPLX( PF2-FMASS**2,MAX(SIGN(FMASS*FWIDTH,PF2),RXZERO))
      SL1= (VC(1)+       VC(4))*FI(1)
     &    +(VC(2)-CXIMAG*VC(3))*FI(2)
      SL2= (VC(2)+CXIMAG*VC(3))*FI(1)
     &    +(VC(1)-       VC(4))*FI(2)
C
      IF ( G(2) .NE. RXZERO ) THEN
         SR1= (VC(1)-       VC(4))*FI(3)
     &       -(VC(2)-CXIMAG*VC(3))*FI(4)
         SR2=-(VC(2)+CXIMAG*VC(3))*FI(3)
     &       +(VC(1)+       VC(4))*FI(4)
C
         FVI(1) = ( G(1)*((PF(0)-PF(3))*SL1 -DCONJG(FVI(6))*SL2)
     &             +G(2)*FMASS*SR1)*D
         FVI(2) = ( G(1)*(      -FVI(6)*SL1 +(PF(0)+PF(3))*SL2)
     &             +G(2)*FMASS*SR2)*D
         FVI(3) = ( G(2)*((PF(0)+PF(3))*SR1 +DCONJG(FVI(6))*SR2)
     &             +G(1)*FMASS*SL1)*D
         FVI(4) = ( G(2)*(       FVI(6)*SR1 +(PF(0)-PF(3))*SR2)
     &             +G(1)*FMASS*SL2)*D
C
      ELSE          
         FVI(1) = G(1)*((PF(0)-PF(3))*SL1 -DCONJG(FVI(6))*SL2)*D
         FVI(2) = G(1)*(      -FVI(6)*SL1 +(PF(0)+PF(3))*SL2)*D
         FVI(3) = G(1)*FMASS*SL1*D
         FVI(4) = G(1)*FMASS*SL2*D
      END IF
C
      RETURN          
      END
C
C ----------------------------------------------------------------------
C
      SUBROUTINE FVOXXX(FO,VC,G,FMASS,FWIDTH , FVO)
C
C this subroutine computes an off-shell fermion wavefunction from a     
C flowing-out external fermion and a vector boson.                      
C                                                                       
C input:                                                                
C       complex fo(6)          : flow-out fermion                   <fo|
C       complex vc(6)          : input    vector                      v 
C       real    g(2)           : coupling constants                  gvf
C       real    fmass          : mass  of output fermion f'             
C       real    fwidth         : width of output fermion f'             
C                                                                       
C output:                                                               
C       complex fvo(6)         : off-shell fermion             <fo,v,f'|
C
#if defined(CERNLIB_IMPNONE)
      IMPLICIT NONE
#endif
      COMPLEX*16 FO(6),VC(6),FVO(6),SL1,SL2,SR1,SR2,D
      REAL*8    G(2),PF(0:3),FMASS,FWIDTH,PF2
      REAL*8 RXZERO, RXONE
      PARAMETER( RXZERO=0.0D0, RXONE=1.0D0 )
      COMPLEX*16 CXIMAG
      LOGICAL FIRST
      SAVE CXIMAG,FIRST
      DATA FIRST/.TRUE./
C
C          Fix compilation with g77
      IF(FIRST) THEN
        FIRST=.FALSE.
        CXIMAG=DCMPLX( RXZERO, RXONE )
      ENDIF
C
      FVO(5) = FO(5)+VC(5)
      FVO(6) = FO(6)+VC(6)
C
      PF(0)=DBLE( FVO(5))
      PF(1)=DBLE( FVO(6))
      PF(2)=DIMAG(FVO(6))
      PF(3)=DIMAG(FVO(5))
      PF2=PF(0)**2-(PF(1)**2+PF(2)**2+PF(3)**2)
C
      D=-RXONE/DCMPLX( PF2-FMASS**2,MAX(SIGN(FMASS*FWIDTH,PF2),RXZERO))
      SL1= (VC(1)+       VC(4))*FO(3)
     &    +(VC(2)+CXIMAG*VC(3))*FO(4)
      SL2= (VC(2)-CXIMAG*VC(3))*FO(3)
     &    +(VC(1)-       VC(4))*FO(4)
C
      IF ( G(2) .NE. RXZERO ) THEN
         SR1= (VC(1)-       VC(4))*FO(1)
     &       -(VC(2)+CXIMAG*VC(3))*FO(2)
         SR2=-(VC(2)-CXIMAG*VC(3))*FO(1)
     &       +(VC(1)+       VC(4))*FO(2)
C
         FVO(1) = ( G(2)*( (PF(0)+PF(3))*SR1        +FVO(6)*SR2)
     &             +G(1)*FMASS*SL1)*D
         FVO(2) = ( G(2)*( DCONJG(FVO(6))*SR1 +(PF(0)-PF(3))*SR2)
     &             +G(1)*FMASS*SL2)*D
         FVO(3) = ( G(1)*( (PF(0)-PF(3))*SL1        -FVO(6)*SL2)
     &             +G(2)*FMASS*SR1)*D
         FVO(4) = ( G(1)*(-DCONJG(FVO(6))*SL1 +(PF(0)+PF(3))*SL2)
     &             +G(2)*FMASS*SR2)*D
C
      ELSE
         FVO(1) = G(1)*FMASS*SL1*D
         FVO(2) = G(1)*FMASS*SL2*D
         FVO(3) = G(1)*( (PF(0)-PF(3))*SL1        -FVO(6)*SL2)*D
         FVO(4) = G(1)*(-DCONJG(FVO(6))*SL1 +(PF(0)+PF(3))*SL2)*D
      END IF
C
      RETURN          
      END
C
C ----------------------------------------------------------------------
C
      SUBROUTINE GGGGXX(WM,W31,WP,W32,G, VERTEX)
C
C this subroutine computes an amplitude of the four-point coupling of   
C the w-, w+ and two w3/z/a.  the amplitude includes the contributions  
C of w exchange diagrams.  the internal w propagator is given in unitary
C gauge.  if one sets wmass=0.0, then the gggg vertex is given (see sect
C 2.9.1 of the manual).
C                                                                       
C input:                                                                
C       complex wm(0:3)        : flow-out w-                         wm 
C       complex w31(0:3)       : first    w3/z/a                     w31
C       complex wp(0:3)        : flow-out w+                         wp 
C       complex w32(0:3)       : second   w3/z/a                     w32
C       real    g              : coupling of w31 with w-/w+             
C                                                  (see the table below)
C                                                                       
C the possible sets of the inputs are as follows:                       
C   -------------------------------------------                         
C   |  wm  |  w31 |  wp  |  w32 |  g31 |  g32 |                         
C   -------------------------------------------                         
C   |  w-  |  w3  |  w+  |  w3  |  gw  |  gw  |                         
C   |  w-  |  w3  |  w+  |  z   |  gw  | gwwz |                         
C   |  w-  |  w3  |  w+  |  a   |  gw  | gwwa |                         
C   |  w-  |  z   |  w+  |  z   | gwwz | gwwz |                         
C   |  w-  |  z   |  w+  |  a   | gwwz | gwwa |                         
C   |  w-  |  a   |  w+  |  a   | gwwa | gwwa |                         
C   -------------------------------------------                         
C where all the bosons are defined by the flowing-out quantum number.   
C                                                                       
C output:                                                               
C       complex vertex         : amplitude          gamma(wm,w31,wp,w32)
C
#if defined(CERNLIB_IMPNONE)
      IMPLICIT NONE
#endif
      COMPLEX*16    WM(6),W31(6),WP(6),W32(6),VERTEX
      COMPLEX*16 DV1(0:3),DV2(0:3),DV3(0:3),DV4(0:3),
     &           DVERTX,V12,V13,V14,V23,V24,V34
      REAL*8       PWM(0:3),PW31(0:3),PWP(0:3),PW32(0:3),G
      REAL*8     DP1(0:3),DP2(0:3),DP3(0:3),DP4(0:3)
      REAL*8 RXZERO, RXONE
      PARAMETER( RXZERO=0.0D0, RXONE=1.0D0 )
C
      PWM(0)=DBLE( WM(5))
      PWM(1)=DBLE( WM(6))
      PWM(2)=DIMAG(WM(6))
      PWM(3)=DIMAG(WM(5))
      PWP(0)=DBLE( WP(5))
      PWP(1)=DBLE( WP(6))
      PWP(2)=DIMAG(WP(6))
      PWP(3)=DIMAG(WP(5))
      PW31(0)=DBLE( W31(5))
      PW31(1)=DBLE( W31(6))
      PW31(2)=DIMAG(W31(6))
      PW31(3)=DIMAG(W31(5))
      PW32(0)=DBLE( W32(5))
      PW32(1)=DBLE( W32(6))
      PW32(2)=DIMAG(W32(6))
      PW32(3)=DIMAG(W32(5))
C
      DV1(0)=DCMPLX(WM(1))
      DV1(1)=DCMPLX(WM(2))
      DV1(2)=DCMPLX(WM(3))
      DV1(3)=DCMPLX(WM(4))
      DP1(0)=DBLE(PWM(0))
      DP1(1)=DBLE(PWM(1))
      DP1(2)=DBLE(PWM(2))
      DP1(3)=DBLE(PWM(3))
      DV2(0)=DCMPLX(W31(1))
      DV2(1)=DCMPLX(W31(2))
      DV2(2)=DCMPLX(W31(3))
      DV2(3)=DCMPLX(W31(4))
      DP2(0)=DBLE(PW31(0))
      DP2(1)=DBLE(PW31(1))
      DP2(2)=DBLE(PW31(2))
      DP2(3)=DBLE(PW31(3))
      DV3(0)=DCMPLX(WP(1))
      DV3(1)=DCMPLX(WP(2))
      DV3(2)=DCMPLX(WP(3))
      DV3(3)=DCMPLX(WP(4))
      DP3(0)=DBLE(PWP(0))
      DP3(1)=DBLE(PWP(1))
      DP3(2)=DBLE(PWP(2))
      DP3(3)=DBLE(PWP(3))
      DV4(0)=DCMPLX(W32(1))
      DV4(1)=DCMPLX(W32(2))
      DV4(2)=DCMPLX(W32(3))
      DV4(3)=DCMPLX(W32(4))
      DP4(0)=DBLE(PW32(0))
      DP4(1)=DBLE(PW32(1))
      DP4(2)=DBLE(PW32(2))
      DP4(3)=DBLE(PW32(3))
C
      V12= DV1(0)*DV2(0)-DV1(1)*DV2(1)-DV1(2)*DV2(2)-DV1(3)*DV2(3)
      V13= DV1(0)*DV3(0)-DV1(1)*DV3(1)-DV1(2)*DV3(2)-DV1(3)*DV3(3)
      V14= DV1(0)*DV4(0)-DV1(1)*DV4(1)-DV1(2)*DV4(2)-DV1(3)*DV4(3)
      V23= DV2(0)*DV3(0)-DV2(1)*DV3(1)-DV2(2)*DV3(2)-DV2(3)*DV3(3)
      V24= DV2(0)*DV4(0)-DV2(1)*DV4(1)-DV2(2)*DV4(2)-DV2(3)*DV4(3)
      V34= DV3(0)*DV4(0)-DV3(1)*DV4(1)-DV3(2)*DV4(2)-DV3(3)*DV4(3)

         DVERTX = V14*V23 -V13*V24
C
         VERTEX = DCMPLX( DVERTX ) * (G*G)
C
      RETURN
      END
C
C ======================================================================
C
      SUBROUTINE GGGXXX(WM,WP,W3,G , VERTEX)
C
C this subroutine computes an amplitude of the three-point coupling of  
C the gauge bosons.                                                     
C                                                                       
C input:                                                                
C       complex wm(6)          : vector               flow-out w-       
C       complex wp(6)          : vector               flow-out w+       
C       complex w3(6)          : vector               j3 or a    or z   
C       real    g              : coupling constant    gw or gwwa or gwwz
C                                                                       
C output:                                                               
C       complex vertex         : amplitude               gamma(wm,wp,w3)
C
#if defined(CERNLIB_IMPNONE)
      IMPLICIT NONE
#endif
      COMPLEX*16 WM(6),WP(6),W3(6),VERTEX, 
     &        XV1,XV2,XV3,V12,V23,V31,P12,P13,P21,P23,P31,P32
      REAL*8    PWM(0:3),PWP(0:3),PW3(0:3),G
      REAL*8 RXZERO, RTENTH
      PARAMETER( RXZERO=0.0D0, RTENTH=0.1D0 )
C
      PWM(0)=DBLE( WM(5))
      PWM(1)=DBLE( WM(6))
      PWM(2)=DIMAG(WM(6))
      PWM(3)=DIMAG(WM(5))
      PWP(0)=DBLE( WP(5))
      PWP(1)=DBLE( WP(6))
      PWP(2)=DIMAG(WP(6))
      PWP(3)=DIMAG(WP(5))
      PW3(0)=DBLE( W3(5))
      PW3(1)=DBLE( W3(6))
      PW3(2)=DIMAG(W3(6))
      PW3(3)=DIMAG(W3(5))
C
      V12=WM(1)*WP(1)-WM(2)*WP(2)-WM(3)*WP(3)-WM(4)*WP(4)
      V23=WP(1)*W3(1)-WP(2)*W3(2)-WP(3)*W3(3)-WP(4)*W3(4)
      V31=W3(1)*WM(1)-W3(2)*WM(2)-W3(3)*WM(3)-W3(4)*WM(4)
      XV1=RXZERO
      XV2=RXZERO
      XV3=RXZERO
      IF ( ABS(WM(1)) .NE. RXZERO ) THEN
         IF (ABS(WM(1)).GE.MAX(ABS(WM(2)),ABS(WM(3)),ABS(WM(4)))
     $        *RTENTH)
     &      XV1=PWM(0)/WM(1)
      ENDIF
      IF ( ABS(WP(1)) .NE. RXZERO) THEN
         IF (ABS(WP(1)).GE.MAX(ABS(WP(2)),ABS(WP(3)),ABS(WP(4)))
     $        *RTENTH)
     &      XV2=PWP(0)/WP(1)
      ENDIF
      IF ( ABS(W3(1)) .NE. RXZERO) THEN
         IF ( ABS(W3(1)).GE.MAX(ABS(W3(2)),ABS(W3(3)),ABS(W3(4)))
     $        *RTENTH)
     &      XV3=PW3(0)/W3(1)
      ENDIF
      P12= (PWM(0)-XV1*WM(1))*WP(1)-(PWM(1)-XV1*WM(2))*WP(2)
     &    -(PWM(2)-XV1*WM(3))*WP(3)-(PWM(3)-XV1*WM(4))*WP(4)
      P13= (PWM(0)-XV1*WM(1))*W3(1)-(PWM(1)-XV1*WM(2))*W3(2)
     &    -(PWM(2)-XV1*WM(3))*W3(3)-(PWM(3)-XV1*WM(4))*W3(4)
      P21= (PWP(0)-XV2*WP(1))*WM(1)-(PWP(1)-XV2*WP(2))*WM(2)
     &    -(PWP(2)-XV2*WP(3))*WM(3)-(PWP(3)-XV2*WP(4))*WM(4)
      P23= (PWP(0)-XV2*WP(1))*W3(1)-(PWP(1)-XV2*WP(2))*W3(2)
     &    -(PWP(2)-XV2*WP(3))*W3(3)-(PWP(3)-XV2*WP(4))*W3(4)
      P31= (PW3(0)-XV3*W3(1))*WM(1)-(PW3(1)-XV3*W3(2))*WM(2)
     &    -(PW3(2)-XV3*W3(3))*WM(3)-(PW3(3)-XV3*W3(4))*WM(4)
      P32= (PW3(0)-XV3*W3(1))*WP(1)-(PW3(1)-XV3*W3(2))*WP(2)
     &    -(PW3(2)-XV3*W3(3))*WP(3)-(PW3(3)-XV3*W3(4))*WP(4)
C
      VERTEX = -(V12*(P13-P23)+V23*(P21-P31)+V31*(P32-P12))*G
C
      RETURN
      END
      SUBROUTINE HIOXXX(FI,FO,GC,SMASS,SWIDTH , HIO)
C
C this subroutine computes an off-shell scalar current from an external
C fermion pair.
C       
C input:
C       complex fi(6)          : flow-in  fermion                   |fi>
C       complex fo(6)          : flow-out fermion                   <fo|
C       complex gc(2)          : coupling constants                 gchf
C       real    smass          : mass  of output scalar s
C       real    swidth         : width of output scalar s
C       
C output:
C       complex hio(3)         : scalar current             j(<fi|s|fo>)
C
#if defined(CERNLIB_IMPNONE)
      IMPLICIT NONE
#endif
      COMPLEX*16 FI(6),FO(6),HIO(3),GC(2),DN
      REAL*8  Q(0:3),SMASS,SWIDTH,Q2
C       
      HIO(2) = FO(5)-FI(5)
      HIO(3) = FO(6)-FI(6)
C       
      Q(0)=DBLE( HIO(2))
      Q(1)=DBLE( HIO(3))
      Q(2)=DIMAG(HIO(3))
      Q(3)=DIMAG(HIO(2))
      Q2=Q(0)**2-(Q(1)**2+Q(2)**2+Q(3)**2)
C
      DN=-DCMPLX(Q2-SMASS**2,DMAX1(DSIGN(SMASS*SWIDTH,Q2),0.D0))
C
      HIO(1) = ( GC(1)*(FO(1)*FI(1)+FO(2)*FI(2))
     &          +GC(2)*(FO(3)*FI(3)+FO(4)*FI(4)) )/DN
C
      RETURN
      END
C ----------------------------------------------------------------------
C
      SUBROUTINE HSSSXX(S1,S2,S3,G,SMASS,SWIDTH , HSSS)
C
C This subroutine computes an off-shell scalar current from the four-   
C scalar coupling.                                                      
C                                                                       
C INPUT:                                                                
C       complex S1(3)          : first  scalar                        S1
C       complex S2(3)          : second scalar                        S2
C       complex S3(3)          : third  scalar                        S3
C       real    G              : coupling constant                 GHHHH
C       real    SMASS          : mass  of OUTPUT scalar S'              
C       real    SWIDTH         : width of OUTPUT scalar S'              
C                                                                       
C OUTPUT:                                                               
C       complex HSSS(3)        : scalar current           J(S':S1,S2,S3)
C
#if defined(CERNLIB_IMPNONE)
      IMPLICIT NONE
#endif
      COMPLEX*16 S1(3),S2(3),S3(3),HSSS(3),DG
      REAL*8     Q(0:3),G,SMASS,SWIDTH,Q2
C
      HSSS(2) = S1(2)+S2(2)+S3(2)
      HSSS(3) = S1(3)+S2(3)+S3(3)
C
      Q(0)=DBLE( HSSS(2))
      Q(1)=DBLE( HSSS(3))
      Q(2)=DIMAG(HSSS(3))
      Q(3)=DIMAG(HSSS(2))
      Q2=Q(0)**2-(Q(1)**2+Q(2)**2+Q(3)**2)
C
      DG=-G/DCMPLX( Q2-SMASS**2,MAX(SIGN(SMASS*SWIDTH ,Q2),0.D0))
C
      HSSS(1) = DG * S1(1)*S2(1)*S3(1)
C
      RETURN
      END
C ----------------------------------------------------------------------
C
      SUBROUTINE HSSXXX(S1,S2,G,SMASS,SWIDTH , HSS)
C
C This subroutine computes an off-shell scalar current from the three-  
C scalar coupling.                                                      
C                                                                       
C INPUT:                                                                
C       complex S1(3)          : first  scalar                        S1
C       complex S2(3)          : second scalar                        S2
C       real    G              : coupling constant                  GHHH
C       real    SMASS          : mass  of OUTPUT scalar S'              
C       real    SWIDTH         : width of OUTPUT scalar S'              
C                                                                       
C OUTPUT:                                                               
C       complex HSS(3)         : scalar current              J(S':S1,S2)
C
#if defined(CERNLIB_IMPNONE)
      IMPLICIT NONE
#endif
      COMPLEX*16 S1(3),S2(3),HSS(3),DG
      REAL*8  Q(0:3),G,SMASS,SWIDTH,Q2
C
      HSS(2) = S1(2)+S2(2)
      HSS(3) = S1(3)+S2(3)
C
      Q(0)=DBLE( HSS(2))
      Q(1)=DBLE( HSS(3))
      Q(2)=DIMAG(HSS(3))
      Q(3)=DIMAG(HSS(2))
      Q2=Q(0)**2-(Q(1)**2+Q(2)**2+Q(3)**2)
C
      DG=-G/DCMPLX( Q2-SMASS**2, MAX(SIGN(SMASS*SWIDTH ,Q2),0.D0))
C
      HSS(1) = DG*S1(1)*S2(1)
C
      RETURN
      END
C
C ======================================================================
C ----------------------------------------------------------------------
C
      SUBROUTINE HVSXXX(VC,SC,G,SMASS,SWIDTH , HVS)
C
C this subroutine computes an off-shell scalar current from the vector- 
C scalar-scalar coupling.  the coupling is absent in the minimal sm in  
C unitary gauge.                                                        
C                                                                       
C input:                                                                
C       complex vc(6)          : input vector                          v
C       complex sc(3)          : input scalar                          s
C       complex g              : coupling constant (s charge)           
C       real    smass          : mass  of output scalar s'              
C       real    swidth         : width of output scalar s'              
C                                                                       
C examples of the coupling constant g for susy particles are as follows:
C   -----------------------------------------------------------         
C   |    s1    | (q,i3) of s1  ||   v=a   |   v=z   |   v=w   |         
C   -----------------------------------------------------------         
C   | nu~_l    | (  0  , +1/2) ||   ---   |  gzn(1) |  gwf(1) |         
C   | e~_l     | ( -1  , -1/2) ||  gal(1) |  gzl(1) |  gwf(1) |         
C   | u~_l     | (+2/3 , +1/2) ||  gau(1) |  gzu(1) |  gwf(1) |         
C   | d~_l     | (-1/3 , -1/2) ||  gad(1) |  gzd(1) |  gwf(1) |         
C   -----------------------------------------------------------         
C   | e~_r-bar | ( +1  ,  0  ) || -gal(2) | -gzl(2) | -gwf(2) |         
C   | u~_r-bar | (-2/3 ,  0  ) || -gau(2) | -gzu(2) | -gwf(2) |         
C   | d~_r-bar | (+1/3 ,  0  ) || -gad(2) | -gzd(2) | -gwf(2) |         
C   -----------------------------------------------------------         
C where the sc charge is defined by the flowing-out quantum number.     
C                                                                       
C output:                                                               
C       complex hvs(3)         : scalar current                j(s':v,s)
C
#if defined(CERNLIB_IMPNONE)
      IMPLICIT NONE
#endif
      COMPLEX*16 VC(6),SC(3),HVS(3),DG,QVV,QPV,G
      REAL*8    QV(0:3),QP(0:3),QA(0:3),SMASS,SWIDTH,Q2
C
      HVS(2) = VC(5)+SC(2)
      HVS(3) = VC(6)+SC(3)
C
      QV(0)=DBLE(  VC(5))
      QV(1)=DBLE(  VC(6))
      QV(2)=DIMAG( VC(6))
      QV(3)=DIMAG( VC(5))
      QP(0)=DBLE(  SC(2))
      QP(1)=DBLE(  SC(3))
      QP(2)=DIMAG( SC(3))
      QP(3)=DIMAG( SC(2))
      QA(0)=DBLE( HVS(2))
      QA(1)=DBLE( HVS(3))
      QA(2)=DIMAG(HVS(3))
      QA(3)=DIMAG(HVS(2))
      Q2=QA(0)**2-(QA(1)**2+QA(2)**2+QA(3)**2)
C
      DG=-G/DCMPLX( Q2-SMASS**2 , MAX(DSIGN( SMASS*SWIDTH ,Q2),0D0) )
      QVV=QV(0)*VC(1)-QV(1)*VC(2)-QV(2)*VC(3)-QV(3)*VC(4)
      QPV=QP(0)*VC(1)-QP(1)*VC(2)-QP(2)*VC(3)-QP(3)*VC(4)
C
      HVS(1) = DG*(2D0*QPV+QVV)*SC(1)
C
      RETURN
      END
C
C ----------------------------------------------------------------------
C
      SUBROUTINE HVVXXX(V1,V2,G,SMASS,SWIDTH , HVV)
C
C this subroutine computes an off-shell scalar current from the vector- 
C vector-scalar coupling.                                               
C                                                                       
C input:                                                                
C       complex v1(6)          : first  vector                        v1
C       complex v2(6)          : second vector                        v2
C       real    g              : coupling constant                  gvvh
C       real    smass          : mass  of output scalar s               
C       real    swidth         : width of output scalar s               
C                                                                       
C output:                                                               
C       complex hvv(3)         : off-shell scalar current     j(s:v1,v2)
C
#if defined(CERNLIB_IMPNONE)
      IMPLICIT NONE
#endif
      COMPLEX*16 V1(6),V2(6),HVV(3),DG
      REAL*8    Q(0:3),G,SMASS,SWIDTH,Q2
      REAL*8 RXZERO
      PARAMETER( RXZERO=0.0D0 )
C
      HVV(2) = V1(5)+V2(5)
      HVV(3) = V1(6)+V2(6)
C
      Q(0)=DBLE( HVV(2))
      Q(1)=DBLE( HVV(3))
      Q(2)=DIMAG(HVV(3))
      Q(3)=DIMAG(HVV(2))
      Q2=Q(0)**2-(Q(1)**2+Q(2)**2+Q(3)**2)
C
      DG=-G/DCMPLX( Q2-SMASS**2 , MAX(SIGN( SMASS*SWIDTH ,Q2),RXZERO) )
C
      HVV(1) = DG*(V1(1)*V2(1)-V1(2)*V2(2)-V1(3)*V2(3)-V1(4)*V2(4))
C
      RETURN
      END
C
C ======================================================================
C
      SUBROUTINE IOSXXX(FI,FO,SC,GC , VERTEX)
C
C This subroutine computes an amplitude of the fermion-fermion-scalar   
C coupling.                                                             
C                                                                       
C INPUT:                                                                
C       complex FI(6)          : flow-in  fermion                   |FI>
C       complex FO(6)          : flow-out fermion                   <FO|
C       complex SC(3)          : input    scalar                      S 
C       complex GC(2)          : coupling constants                 GCHF
C                                                                       
C OUTPUT:                                                               
C       complex VERTEX         : amplitude                     <FO|S|FI>
C
#if defined(CERNLIB_IMPNONE)
      IMPLICIT NONE
#endif
      COMPLEX*16 FI(6),FO(6),SC(3),GC(2),VERTEX
C
      VERTEX = SC(1)*( GC(1)*(FI(1)*FO(1)+FI(2)*FO(2))
     &                +GC(2)*(FI(3)*FO(3)+FI(4)*FO(4)) )
C
      RETURN          
      END
C
C ======================================================================
C
      SUBROUTINE IOVXXX(FI,FO,VC,G , VERTEX)
C
C this subroutine computes an amplitude of the fermion-fermion-vector   
C coupling.                                                             
C                                                                       
C input:                                                                
C       complex fi(6)          : flow-in  fermion                   |fi>
C       complex fo(6)          : flow-out fermion                   <fo|
C       complex vc(6)          : input    vector                      v 
C       real    g(2)           : coupling constants                  gvf
C                                                                       
C output:                                                               
C       complex vertex         : amplitude                     <fo|v|fi>
C
#if defined(CERNLIB_IMPNONE)
      IMPLICIT NONE
#endif
      COMPLEX*16 FI(6),FO(6),VC(6),VERTEX
      REAL*8    G(2)
      REAL*8 RXZERO, RXONE
      PARAMETER( RXZERO=0.0D0, RXONE=1.0D0 )
      COMPLEX*16 CXIMAG
      LOGICAL FIRST
      SAVE CXIMAG,FIRST
      DATA FIRST/.TRUE./
C
C          Fix compilation with g77
      IF(FIRST) THEN
        FIRST=.FALSE.
        CXIMAG=DCMPLX( RXZERO, RXONE )
      ENDIF
C

      VERTEX =  G(1)*( (FO(3)*FI(1)+FO(4)*FI(2))*VC(1)
     &                +(FO(3)*FI(2)+FO(4)*FI(1))*VC(2)
     &                -(FO(3)*FI(2)-FO(4)*FI(1))*VC(3)*CXIMAG
     &                +(FO(3)*FI(1)-FO(4)*FI(2))*VC(4)        )
C
      IF ( G(2) .NE. RXZERO ) THEN
         VERTEX = VERTEX
     &          + G(2)*( (FO(1)*FI(3)+FO(2)*FI(4))*VC(1)
     &                  -(FO(1)*FI(4)+FO(2)*FI(3))*VC(2)
     &                  +(FO(1)*FI(4)-FO(2)*FI(3))*VC(3)*CXIMAG
     &                  -(FO(1)*FI(3)-FO(2)*FI(4))*VC(4)        )
      END IF
C
      RETURN
      END
C
C	Subroutine returns the desired fermion or
C	anti-fermion spinor. ie., |f>
C	A replacement for the HELAS routine IXXXXX
C
C	Adam Duff,  1992 August 31
C	<duff@phenom.physics.wisc.edu>
C
      SUBROUTINE IXXXXX(P,FMASS,NHEL,NSF,FI)
C          P        IN: FOUR VECTOR MOMENTUM
C          FMASS    IN: FERMION MASS
C          NHEL     IN: SPINOR HELICITY, -1 OR 1
C          NSF      IN: -1=ANTIFERMION, 1=FERMION
C          FI       OUT: FERMION WAVEFUNCTION
C
C declare input/output variables
C
#if defined(CERNLIB_IMPNONE)
      IMPLICIT NONE
#endif
      COMPLEX*16 FI(6)
      INTEGER*4 NHEL, NSF
      REAL*8 P(0:3), FMASS
      REAL*8 RXZERO, RXONE, RXTWO
      PARAMETER( RXZERO=0.0D0, RXONE=1.0D0, RXTWO=2.0D0 )
      REAL*8 PLAT, PABS, OMEGAP, OMEGAM, RS2PA, SPAZ
      COMPLEX*16 CXZERO
C
C declare local variables
C
      LOGICAL FIRST
      SAVE CXZERO,FIRST
      DATA FIRST/.TRUE./
C
C          Fix compilation with g77
      IF(FIRST) THEN
        FIRST=.FALSE.
        CXZERO=DCMPLX( RXZERO, RXZERO )
      ENDIF
C
C define kinematic parameters
C
      FI(5) = DCMPLX( P(0), P(3) ) * NSF
      FI(6) = DCMPLX( P(1), P(2) ) * NSF
      PLAT = SQRT( P(1)**2 + P(2)**2 )
      PABS = SQRT( P(1)**2 + P(2)**2 + P(3)**2 )
      OMEGAP = SQRT( P(0) + PABS )
C
C do massive fermion case
C
      IF ( FMASS .NE. RXZERO ) THEN
         OMEGAM = FMASS / OMEGAP
         IF ( NSF .EQ. 1 ) THEN
            IF ( NHEL .EQ. 1 ) THEN
               IF ( P(3) .GE. RXZERO ) THEN
                  IF ( PLAT .EQ. RXZERO ) THEN
                     FI(1) = DCMPLX( OMEGAM, RXZERO )
                     FI(2) = CXZERO
                     FI(3) = DCMPLX( OMEGAP, RXZERO )
                     FI(4) = CXZERO
                  ELSE
                     RS2PA = RXONE / SQRT( RXTWO * PABS )
                     SPAZ = SQRT( PABS + P(3) )
                     FI(1) = OMEGAM * RS2PA
     &                     * DCMPLX( SPAZ, RXZERO )
                     FI(2) = OMEGAM * RS2PA / SPAZ
     &                     * DCMPLX( P(1), P(2) )
                     FI(3) = OMEGAP * RS2PA
     &                     * DCMPLX( SPAZ, RXZERO )
                     FI(4) = OMEGAP * RS2PA / SPAZ
     &                     * DCMPLX( P(1), P(2) )
                  END IF
               ELSE
                  IF ( PLAT .EQ. RXZERO ) THEN
                     FI(1) = CXZERO
                     FI(2) = DCMPLX( OMEGAM, RXZERO )
                     FI(3) = CXZERO
                     FI(4) = DCMPLX( OMEGAP, RXZERO )
                  ELSE
                     RS2PA = RXONE / SQRT( RXTWO * PABS )
                     SPAZ = SQRT( PABS - P(3) )
                     FI(1) = OMEGAM * RS2PA / SPAZ
     &                     * DCMPLX( PLAT, RXZERO )
                     FI(2) = OMEGAM * RS2PA * SPAZ / PLAT
     &                     * DCMPLX( P(1), P(2) )
                     FI(3) = OMEGAP * RS2PA / SPAZ
     &                     * DCMPLX( PLAT, RXZERO )
                     FI(4) = OMEGAP * RS2PA * SPAZ / PLAT
     &                     * DCMPLX( P(1), P(2) )
                  END IF
               END IF
            ELSE IF ( NHEL .EQ. -1 ) THEN
               IF ( P(3) .GE. RXZERO ) THEN
                  IF ( PLAT .EQ. RXZERO ) THEN
                     FI(1) = CXZERO
                     FI(2) = DCMPLX( OMEGAP, RXZERO )
                     FI(3) = CXZERO
                     FI(4) = DCMPLX( OMEGAM, RXZERO )
                  ELSE
                     RS2PA = RXONE / SQRT( RXTWO * PABS )
                     SPAZ = SQRT( PABS + P(3) )
                     FI(1) = OMEGAP * RS2PA / SPAZ
     &                     * DCMPLX( -P(1), P(2) )
                     FI(2) = OMEGAP * RS2PA
     &                     * DCMPLX( SPAZ, RXZERO )
                     FI(3) = OMEGAM * RS2PA / SPAZ
     &                     * DCMPLX( -P(1), P(2) )
                     FI(4) = OMEGAM * RS2PA
     &                     * DCMPLX( SPAZ, RXZERO )
                  END IF
               ELSE
                  IF ( PLAT .EQ. RXZERO ) THEN
                     FI(1) = DCMPLX( -OMEGAP, RXZERO )
                     FI(2) = CXZERO
                     FI(3) = DCMPLX( -OMEGAM, RXZERO )
                     FI(4) = CXZERO
                  ELSE
                     RS2PA = RXONE / SQRT( RXTWO * PABS )
                     SPAZ = SQRT( PABS - P(3) )
                     FI(1) = OMEGAP * RS2PA * SPAZ / PLAT
     &                     * DCMPLX( -P(1), P(2) )
                     FI(2) = OMEGAP * RS2PA / SPAZ
     &                     * DCMPLX( PLAT, RXZERO )
                     FI(3) = OMEGAM * RS2PA * SPAZ / PLAT
     &                     * DCMPLX( -P(1), P(2) )
                     FI(4) = OMEGAM * RS2PA / SPAZ
     &                     * DCMPLX( PLAT, RXZERO )
                  END IF
               END IF
            ELSE
               STOP 'IXXXXX:  FERMION HELICITY MUST BE +1,-1'
            END IF
         ELSE IF ( NSF .EQ. -1 ) THEN
            IF ( NHEL .EQ. 1 ) THEN
               IF ( P(3) .GE. RXZERO ) THEN
                  IF ( PLAT .EQ. RXZERO ) THEN
                     FI(1) = CXZERO
                     FI(2) = DCMPLX( -OMEGAP, RXZERO )
                     FI(3) = CXZERO
                     FI(4) = DCMPLX( OMEGAM, RXZERO )
                  ELSE
                     RS2PA = RXONE / SQRT( RXTWO * PABS )
                     SPAZ = SQRT( PABS + P(3) )
                     FI(1) = -OMEGAP * RS2PA / SPAZ
     &                     * DCMPLX( -P(1), P(2) )
                     FI(2) = -OMEGAP * RS2PA
     &                     * DCMPLX( SPAZ, RXZERO )
                     FI(3) = OMEGAM * RS2PA / SPAZ
     &                     * DCMPLX( -P(1), P(2) )
                     FI(4) = OMEGAM * RS2PA
     &                     * DCMPLX( SPAZ, RXZERO )
                  END IF
               ELSE
                  IF ( PLAT .EQ. RXZERO ) THEN
                     FI(1) = DCMPLX( OMEGAP, RXZERO )
                     FI(2) = CXZERO
                     FI(3) = DCMPLX( -OMEGAM, RXZERO )
                     FI(4) = CXZERO
                  ELSE
                     RS2PA = RXONE / SQRT( RXTWO * PABS )
                     SPAZ = SQRT( PABS - P(3) )
                     FI(1) = -OMEGAP * RS2PA * SPAZ / PLAT
     &                     * DCMPLX( -P(1), P(2) )
                     FI(2) = -OMEGAP * RS2PA / SPAZ
     &                     * DCMPLX( PLAT, RXZERO )
                     FI(3) = OMEGAM * RS2PA * SPAZ / PLAT
     &                     * DCMPLX( -P(1), P(2) )
                     FI(4) = OMEGAM * RS2PA / SPAZ
     &                     * DCMPLX( PLAT, RXZERO )
                  END IF
               END IF
            ELSE IF ( NHEL .EQ. -1 ) THEN
               IF ( P(3) .GE. RXZERO ) THEN
                  IF ( PLAT .EQ. RXZERO ) THEN
                     FI(1) = DCMPLX( OMEGAM, RXZERO )
                     FI(2) = CXZERO
                     FI(3) = DCMPLX( -OMEGAP, RXZERO )
                     FI(4) = CXZERO
                  ELSE
                     RS2PA = RXONE / SQRT( RXTWO * PABS )
                     SPAZ = SQRT( PABS + P(3) )
                     FI(1) = OMEGAM * RS2PA
     &                     * DCMPLX( SPAZ, RXZERO )
                     FI(2) = OMEGAM * RS2PA / SPAZ
     &                     * DCMPLX( P(1), P(2) )
                     FI(3) = -OMEGAP * RS2PA
     &                     * DCMPLX( SPAZ, RXZERO )
                     FI(4) = -OMEGAP * RS2PA / SPAZ
     &                     * DCMPLX( P(1), P(2) )
                  END IF
               ELSE
                  IF ( PLAT .EQ. RXZERO ) THEN
                     FI(1) = CXZERO
                     FI(2) = DCMPLX( OMEGAM, RXZERO )
                     FI(3) = CXZERO
                     FI(4) = DCMPLX( -OMEGAP, RXZERO )
                  ELSE
                     RS2PA = RXONE / SQRT( RXTWO * PABS )
                     SPAZ = SQRT( PABS - P(3) )
                     FI(1) = OMEGAM * RS2PA / SPAZ
     &                     * DCMPLX( PLAT, RXZERO )
                     FI(2) = OMEGAM * RS2PA * SPAZ / PLAT
     &                     * DCMPLX( P(1), P(2) )
                     FI(3) = -OMEGAP * RS2PA / SPAZ
     &                     * DCMPLX( PLAT, RXZERO )
                     FI(4) = -OMEGAP * RS2PA * SPAZ / PLAT
     &                     * DCMPLX( P(1), P(2) )
                  END IF
               END IF
            ELSE
               STOP 'IXXXXX:  FERMION HELICITY MUST BE +1,-1'
            END IF
         ELSE
            STOP 'IXXXXX:  FERMION TYPE MUST BE +1,-1'
         END IF
C
C do massless fermion case
C
      ELSE
         IF ( NSF .EQ. 1 ) THEN
            IF ( NHEL .EQ. 1 ) THEN
               IF ( P(3) .GE. RXZERO ) THEN
                  IF ( PLAT .EQ. RXZERO ) THEN
                     FI(1) = CXZERO
                     FI(2) = CXZERO
                     FI(3) = DCMPLX( OMEGAP, RXZERO )
                     FI(4) = CXZERO
                  ELSE
                     SPAZ = SQRT( PABS + P(3) )
                     FI(1) = CXZERO
                     FI(2) = CXZERO
                     FI(3) = DCMPLX( SPAZ, RXZERO )
                     FI(4) = RXONE / SPAZ
     &                     * DCMPLX( P(1), P(2) )
                  END IF
               ELSE
                  IF ( PLAT .EQ. RXZERO ) THEN
                     FI(1) = CXZERO
                     FI(2) = CXZERO
                     FI(3) = CXZERO
                     FI(4) = DCMPLX( OMEGAP, RXZERO )
                  ELSE
                     SPAZ = SQRT( PABS - P(3) )
                     FI(1) = CXZERO
                     FI(2) = CXZERO
                     FI(3) = RXONE / SPAZ
     &                     * DCMPLX( PLAT, RXZERO )
                     FI(4) = SPAZ / PLAT
     &                     * DCMPLX( P(1), P(2) )
                  END IF
               END IF
            ELSE IF ( NHEL .EQ. -1 ) THEN
               IF ( P(3) .GE. RXZERO ) THEN
                  IF ( PLAT .EQ. RXZERO ) THEN
                     FI(1) = CXZERO
                     FI(2) = DCMPLX( OMEGAP, RXZERO )
                     FI(3) = CXZERO
                     FI(4) = CXZERO
                  ELSE
                     SPAZ = SQRT( PABS + P(3) )
                     FI(1) = RXONE / SPAZ
     &                     * DCMPLX( -P(1), P(2) )
                     FI(2) = DCMPLX( SPAZ, RXZERO )
                     FI(3) = CXZERO
                     FI(4) = CXZERO
                  END IF
               ELSE
                  IF ( PLAT .EQ. RXZERO ) THEN
                     FI(1) = DCMPLX( -OMEGAP, RXZERO )
                     FI(2) = CXZERO
                     FI(3) = CXZERO
                     FI(4) = CXZERO
                  ELSE
                     SPAZ = SQRT( PABS - P(3) )
                     FI(1) = SPAZ / PLAT
     &                     * DCMPLX( -P(1), P(2) )
                     FI(2) = RXONE / SPAZ
     &                     * DCMPLX( PLAT, RXZERO )
                     FI(3) = CXZERO
                     FI(4) = CXZERO
                  END IF
               END IF
            ELSE
               STOP 'IXXXXX:  FERMION HELICITY MUST BE +1,-1'
            END IF
         ELSE IF ( NSF .EQ. -1 ) THEN
            IF ( NHEL .EQ. 1 ) THEN
               IF ( P(3) .GE. RXZERO ) THEN
                  IF ( PLAT .EQ. RXZERO ) THEN
                     FI(1) = CXZERO
                     FI(2) = DCMPLX( -OMEGAP, RXZERO )
                     FI(3) = CXZERO
                     FI(4) = CXZERO
                  ELSE
                     SPAZ = SQRT( PABS + P(3) )
                     FI(1) = -RXONE / SPAZ
     &                     * DCMPLX( -P(1), P(2) )
                     FI(2) = DCMPLX( -SPAZ, RXZERO )
                     FI(3) = CXZERO
                     FI(4) = CXZERO
                  END IF
               ELSE
                  IF ( PLAT .EQ. RXZERO ) THEN
                     FI(1) = DCMPLX( OMEGAP, RXZERO )
                     FI(2) = CXZERO
                     FI(3) = CXZERO
                     FI(4) = CXZERO
                  ELSE
                     SPAZ = SQRT( PABS - P(3) )
                     FI(1) = -SPAZ / PLAT
     &                     * DCMPLX( -P(1), P(2) )
                     FI(2) = -RXONE / SPAZ
     &                     * DCMPLX( PLAT, RXZERO )
                     FI(3) = CXZERO
                     FI(4) = CXZERO
                  END IF
               END IF
            ELSE IF ( NHEL .EQ. -1 ) THEN
               IF ( P(3) .GE. RXZERO ) THEN
                  IF ( PLAT .EQ. RXZERO ) THEN
                     FI(1) = CXZERO
                     FI(2) = CXZERO
                     FI(3) = DCMPLX( -OMEGAP, RXZERO )
                     FI(4) = CXZERO
                  ELSE
                     SPAZ = SQRT( PABS + P(3) )
                     FI(1) = CXZERO
                     FI(2) = CXZERO
                     FI(3) = DCMPLX( -SPAZ, RXZERO )
                     FI(4) = -RXONE / SPAZ
     &                     * DCMPLX( P(1), P(2) )
                  END IF
               ELSE
                  IF ( PLAT .EQ. RXZERO ) THEN
                     FI(1) = CXZERO
                     FI(2) = CXZERO
                     FI(3) = CXZERO
                     FI(4) = DCMPLX( -OMEGAP, RXZERO )
                  ELSE
                     SPAZ = SQRT( PABS - P(3) )
                     FI(1) = CXZERO
                     FI(2) = CXZERO
                     FI(3) = -RXONE / SPAZ
     &                     * DCMPLX( PLAT, RXZERO )
                     FI(4) = -SPAZ / PLAT
     &                     * DCMPLX( P(1), P(2) )
                  END IF
               END IF
            ELSE
               STOP 'IXXXXX:  FERMION HELICITY MUST BE +1,-1'
            END IF
         ELSE
            STOP 'IXXXXX:  FERMION TYPE MUST BE +1,-1'
         END IF
      END IF
C
C done
C
      RETURN
      END
C
C ----------------------------------------------------------------------
C
      SUBROUTINE J3XXXX(FI,FO,GAF,GZF,ZMASS,ZWIDTH , J3)
C
C this subroutine computes the sum of photon and z currents with the    
C suitable weights ( j(w3) = cos(theta_w) j(z) + sin(theta_w) j(a) ).   
C the output j3 is useful as an input of vvvxxx, jvvxxx or w3w3xx.      
C the photon propagator is given in feynman gauge, and the z propagator 
C is given in unitary gauge.                                            
C                                                                       
C input:                                                                
C       complex fi(6)          : flow-in  fermion                   |fi>
C       complex fo(6)          : flow-out fermion                   <fo|
C       real    gaf(2)         : fi couplings with a                 gaf
C       real    gzf(2)         : fi couplings with z                 gzf
C       real    zmass          : mass  of z                             
C       real    zwidth         : width of z                             
C                                                                       
C output:                                                               
C       complex j3(6)          : w3 current             j^mu(<fo|w3|fi>)
C
#if defined(CERNLIB_IMPNONE)
      IMPLICIT NONE
#endif
      COMPLEX*16 FI(6),FO(6),J3(6),
     &        C0L,C1L,C2L,C3L,CSL,C0R,C1R,C2R,C3R,CSR,DZ,DDIF
      REAL*8    GAF(2),GZF(2),Q(0:3),ZMASS,ZWIDTH,ZM2,ZMW,Q2,DA,WW,
     &        CW,SW,GN,GZ3L,GA3L
C
      REAL*8 RXZERO, RXONE
      PARAMETER( RXZERO=0.0D0, RXONE=1.0D0 )
      COMPLEX*16 CXIMAG
      LOGICAL FIRST
      SAVE CXIMAG,FIRST
      DATA FIRST/.TRUE./
C
C          Fix compilation with g77
      IF(FIRST) THEN
        FIRST=.FALSE.
        CXIMAG=DCMPLX( RXZERO, RXONE )
      ENDIF
C
      J3(5) = FO(5)-FI(5)
      J3(6) = FO(6)-FI(6)
C
      Q(0)=-DBLE( J3(5))
      Q(1)=-DBLE( J3(6))
      Q(2)=-DIMAG(J3(6))
      Q(3)=-DIMAG(J3(5))
      Q2=Q(0)**2-(Q(1)**2+Q(2)**2+Q(3)**2)
      ZM2=ZMASS**2
      ZMW=ZMASS*ZWIDTH
C
      DA=RXONE/Q2
      WW=MAX(DSIGN( ZMW ,Q2),RXZERO)
      DZ=RXONE/DCMPLX( Q2-ZM2 , WW )
      DDIF=DCMPLX( -ZM2 , WW )*DA*DZ
C
C ddif is the difference : ddif=da-dz
C  for the running width, use below instead of the above ww,dz and ddif.
C      ww=max( zwidth*q2/zmass ,r_zero)
C      dz=r_one/dcmplx( q2-zm2 , ww )
C      ddif=dcmplx( -zm2 , ww )*da*dz
C
      CW=RXONE/SQRT(RXONE+(GZF(2)/GAF(2))**2)
      SW=SQRT((RXONE-CW)*(RXONE+CW))
      GN=GAF(2)*SW
      GZ3L=GZF(1)*CW
      GA3L=GAF(1)*SW
      C0L=  FO(3)*FI(1)+FO(4)*FI(2)
      C0R=  FO(1)*FI(3)+FO(2)*FI(4)
      C1L=-(FO(3)*FI(2)+FO(4)*FI(1))
      C1R=  FO(1)*FI(4)+FO(2)*FI(3)
      C2L= (FO(3)*FI(2)-FO(4)*FI(1))*CXIMAG
      C2R=(-FO(1)*FI(4)+FO(2)*FI(3))*CXIMAG
      C3L= -FO(3)*FI(1)+FO(4)*FI(2)
      C3R=  FO(1)*FI(3)-FO(2)*FI(4)
      CSL=(Q(0)*C0L-Q(1)*C1L-Q(2)*C2L-Q(3)*C3L)/ZM2
      CSR=(Q(0)*C0R-Q(1)*C1R-Q(2)*C2R-Q(3)*C3R)/ZM2
C
      J3(1) =  GZ3L*DZ*(C0L-CSL*Q(0))+GA3L*C0L*DA
     &       + GN*(C0R*DDIF-CSR*Q(0)*DZ)
      J3(2) =  GZ3L*DZ*(C1L-CSL*Q(1))+GA3L*C1L*DA
     &       + GN*(C1R*DDIF-CSR*Q(1)*DZ)
      J3(3) =  GZ3L*DZ*(C2L-CSL*Q(2))+GA3L*C2L*DA
     &       + GN*(C2R*DDIF-CSR*Q(2)*DZ)
      J3(4) =  GZ3L*DZ*(C3L-CSL*Q(3))+GA3L*C3L*DA
     &       + GN*(C3R*DDIF-CSR*Q(3)*DZ)
C
      RETURN
      END
C
C ----------------------------------------------------------------------
C
      SUBROUTINE JEEXXX(EB,EF,SHLF,CHLF,PHI,NHB,NHF,NSF , JEE)
C
C This subroutine computes an off-shell photon wavefunction emitted from
C the electron or positron beam, with a special care for the small angle
C region.  The momenta are measured in the laboratory frame, where the  
C e- (e+) beam is along the positive (negative) z axis.                 
C                                                                       
C INPUT:                                                                
C       real    EB             : energy (GeV)    of beam  e-/e+         
C       real    EF             : energy (GeV)    of final e-/e+         
C       real    SHLF           : sin(theta/2)    of final e-/e+         
C       real    CHLF           : cos(theta/2)    of final e-/e+         
C       real    PHI            : azimuthal angle of final e-/e+         
C       integer NHB  = -1 or 1 : helicity        of beam  e-/e+         
C       integer NHF  = -1 or 1 : helicity        of final e-/e+         
C       integer NSF  = -1 or 1 : +1 for electron, -1 for positron       
C                                                                       
C OUTPUT:                                                               
C       complex JEE(6)         : off-shell photon          J^mu(<e|A|e>)
C
#if defined(CERNLIB_IMPNONE)
      IMPLICIT NONE
#endif
      COMPLEX*16 JEE(6),COEFF
      REAL*8  CS(2),EB,EF,SHLF,CHLF,PHI,ME,ALPHA,GAL,HI,SF,SFH,X,ME2,Q2,
     &        RFP,RFM,SNP,CSP,RXC,C,S
      INTEGER NHB,NHF,NSF
C
      ME   =0.51099906D-3
      ALPHA=1./128.
      GAL  =SQRT(ALPHA*4.*3.14159265D0)
C
      HI =NHB
      SF =NSF
      SFH=NHB*NSF
      CS((3+NSF)/2)=SHLF
      CS((3-NSF)/2)=CHLF
C CS(1)=CHLF and CS(2)=SHLF for electron
C CS(1)=SHLF and CS(2)=CHLF for positron
      X=EF/EB
      ME2=ME**2
      Q2=-4.*CS(2)**2*(EF*EB-ME2)
     &   +SF*(1.-X)**2/X*(SHLF+CHLF)*(SHLF-CHLF)*ME2
      RFP=(1+NSF)
      RFM=(1-NSF)
      SNP=SIN(PHI)
      CSP=COS(PHI)
C
      IF (NHB.EQ.NHF) THEN
         RXC=2.*X/(1.-X)*CS(1)**2
         COEFF= GAL*2.*EB*SQRT(X)*CS(2)/Q2
     &         *(DCMPLX( RFP )-RFM*DCMPLX( CSP ,-SNP*HI ))*.5
         JEE(1) =  DCMPLX( 0.D0 )
         JEE(2) =  COEFF*DCMPLX( (1.+RXC)*CSP ,-SFH*SNP )
         JEE(3) =  COEFF*DCMPLX( (1.+RXC)*SNP , SFH*CSP )
         JEE(4) =  COEFF*(-SF*RXC/CS(1)*CS(2))
      ELSE
         COEFF= GAL*ME/Q2/SQRT(X)
     &         *(DCMPLX( RFP )+RFM*DCMPLX( CSP , SNP*HI ))*.5*HI
         JEE(1) = -COEFF*(1.+X)*CS(2)*DCMPLX( CSP , SFH*SNP )
         JEE(2) =  COEFF*(1.-X)*CS(1)
         JEE(3) =  JEE(2)*DCMPLX( 0.D0 , SFH )
         JEE(4) =  JEE(1)*SF*(1.-X)/(1.+X)
      ENDIF
C
      C=(CHLF+SHLF)*(CHLF-SHLF)
      S=2.*CHLF*SHLF
C
      JEE(5) = -EB*DCMPLX( 1.-X , SF-X*C )
      JEE(6) =  EB*X*S*DCMPLX( CSP , SNP )
C
      RETURN          
      END
C
C
C ----------------------------------------------------------------------
C
      SUBROUTINE JGGGXX(W1,W2,W3,G, JW3W)
C
C this subroutine computes an off-shell w+, w-, w3, z or photon current 
C from the four-point gauge boson coupling, including the contributions 
C of w exchange diagrams.  the vector propagator is given in feynman    
C gauge for a photon and in unitary gauge for w and z bosons.  if one   
C sets wmass=0.0, then the ggg-->g current is given (see sect 2.9.1 of 
C the manual).                                                          
C                                                                       
C input:                                                                
C       complex w1(6)          : first  vector                        w1
C       complex w2(6)          : second vector                        w2
C       complex w3(6)          : third  vector                        w3
C       real    g             : first  coupling constant               
C                                                  (see the table below)
C                                                                       
C output:                                                               
C       complex jw3w(6)        : w current             j^mu(w':w1,w2,w3)
C
#if defined(CERNLIB_IMPNONE)
      IMPLICIT NONE
#endif
      COMPLEX*16  W1(6),W2(6),W3(6),JW3W(6)
      COMPLEX*16 DW1(0:3),DW2(0:3),DW3(0:3),
     &           JJ(0:3),DV,W32,W13
      REAL*8     P1(0:3),P2(0:3),P3(0:3),Q(0:3),G,DG2,Q2
C
      REAL*8 RXZERO
      PARAMETER( RXZERO=0.0D0 )
C
      JW3W(5) = W1(5)+W2(5)+W3(5)
      JW3W(6) = W1(6)+W2(6)+W3(6)
C
      DW1(0)=DCMPLX(W1(1))
      DW1(1)=DCMPLX(W1(2))
      DW1(2)=DCMPLX(W1(3))
      DW1(3)=DCMPLX(W1(4))
      DW2(0)=DCMPLX(W2(1))
      DW2(1)=DCMPLX(W2(2))
      DW2(2)=DCMPLX(W2(3))
      DW2(3)=DCMPLX(W2(4))
      DW3(0)=DCMPLX(W3(1))
      DW3(1)=DCMPLX(W3(2))
      DW3(2)=DCMPLX(W3(3))
      DW3(3)=DCMPLX(W3(4))
      P1(0)=DBLE(      W1(5))
      P1(1)=DBLE(      W1(6))
      P1(2)=DBLE(DIMAG(W1(6)))
      P1(3)=DBLE(DIMAG(W1(5)))
      P2(0)=DBLE(      W2(5))
      P2(1)=DBLE(      W2(6))
      P2(2)=DBLE(DIMAG(W2(6)))
      P2(3)=DBLE(DIMAG(W2(5)))
      P3(0)=DBLE(      W3(5))
      P3(1)=DBLE(      W3(6))
      P3(2)=DBLE(DIMAG(W3(6)))
      P3(3)=DBLE(DIMAG(W3(5)))
      Q(0)=-(P1(0)+P2(0)+P3(0))
      Q(1)=-(P1(1)+P2(1)+P3(1))
      Q(2)=-(P1(2)+P2(2)+P3(2))
      Q(3)=-(P1(3)+P2(3)+P3(3))

      Q2 =Q(0)**2 -(Q(1)**2 +Q(2)**2 +Q(3)**2)

      DG2=DBLE(G)*DBLE(G)
C
      DV = 1.0D0/DCMPLX( Q2 )

C  for the running width, use below instead of the above dv.
C      dv = 1.0d0/dcmplx( q2 -mv2 , dmax1(dwv*q2/dmv,0.d0) )
C
      W32=DW3(0)*DW2(0)-DW3(1)*DW2(1)-DW3(2)*DW2(2)-DW3(3)*DW2(3)
C
C     
      W13=DW1(0)*DW3(0)-DW1(1)*DW3(1)-DW1(2)*DW3(2)-DW1(3)*DW3(3)
C     
      JJ(0)=DG2*( DW1(0)*W32 - DW2(0)*W13 )
      JJ(1)=DG2*( DW1(1)*W32 - DW2(1)*W13 )
      JJ(2)=DG2*( DW1(2)*W32 - DW2(2)*W13 )
      JJ(3)=DG2*( DW1(3)*W32 - DW2(3)*W13 )
C     
      JW3W(1) = DCMPLX( JJ(0)*DV )
      JW3W(2) = DCMPLX( JJ(1)*DV )
      JW3W(3) = DCMPLX( JJ(2)*DV )
      JW3W(4) = DCMPLX( JJ(3)*DV )
C
      RETURN
      END
C
C ----------------------------------------------------------------------
C
      SUBROUTINE JGGXXX(V1,V2,G, JVV)
C
C this subroutine computes an off-shell vector current from the three-  
C point gauge boson coupling.  the vector propagator is given in feynman
C gauge for a massless vector and in unitary gauge for a massive vector.
C                                                                       
C input:                                                                
C       complex v1(6)          : first  vector                        v1
C       complex v2(6)          : second vector                        v2
C       real    g              : coupling constant (see the table below)
C                                                                       
C output:                                                               
C       complex jvv(6)         : vector current            j^mu(v:v1,v2)
C
#if defined(CERNLIB_IMPNONE)
      IMPLICIT NONE
#endif
      COMPLEX*16 V1(6),V2(6),JVV(6),J12(0:3),
     &        SV1,SV2,V12
      REAL*8    P1(0:3),P2(0:3),Q(0:3),G,GS,S
C
      REAL*8 RXZERO
      PARAMETER( RXZERO=0.0D0 )
C
      JVV(5) = V1(5)+V2(5)
      JVV(6) = V1(6)+V2(6)
C
      P1(0)=DBLE( V1(5))
      P1(1)=DBLE( V1(6))
      P1(2)=DIMAG(V1(6))
      P1(3)=DIMAG(V1(5))
      P2(0)=DBLE( V2(5))
      P2(1)=DBLE( V2(6))
      P2(2)=DIMAG(V2(6))
      P2(3)=DIMAG(V2(5))
      Q(0)=-DBLE( JVV(5))
      Q(1)=-DBLE( JVV(6))
      Q(2)=-DIMAG(JVV(6))
      Q(3)=-DIMAG(JVV(5))
      S=Q(0)**2-(Q(1)**2+Q(2)**2+Q(3)**2)
C
      V12=V1(1)*V2(1)-V1(2)*V2(2)-V1(3)*V2(3)-V1(4)*V2(4)
      SV1= (P2(0)-Q(0))*V1(1) -(P2(1)-Q(1))*V1(2)
     &    -(P2(2)-Q(2))*V1(3) -(P2(3)-Q(3))*V1(4)
      SV2=-(P1(0)-Q(0))*V2(1) +(P1(1)-Q(1))*V2(2)
     &    +(P1(2)-Q(2))*V2(3) +(P1(3)-Q(3))*V2(4)
      J12(0)=(P1(0)-P2(0))*V12 +SV1*V2(1) +SV2*V1(1)
      J12(1)=(P1(1)-P2(1))*V12 +SV1*V2(2) +SV2*V1(2)
      J12(2)=(P1(2)-P2(2))*V12 +SV1*V2(3) +SV2*V1(3)
      J12(3)=(P1(3)-P2(3))*V12 +SV1*V2(4) +SV2*V1(4)
C
      GS=-G/S
C
      JVV(1) = GS*J12(0)
      JVV(2) = GS*J12(1)
      JVV(3) = GS*J12(2)
      JVV(4) = GS*J12(3)
C
      RETURN
      END
C
C ----------------------------------------------------------------------
C
      SUBROUTINE JIOXXX(FI,FO,G,VMASS,VWIDTH , JIO)
C
C this subroutine computes an off-shell vector current from an external 
C fermion pair.  the vector boson propagator is given in feynman gauge  
C for a massless vector and in unitary gauge for a massive vector.      
C                                                                       
C input:                                                                
C       complex fi(6)          : flow-in  fermion                   |fi>
C       complex fo(6)          : flow-out fermion                   <fo|
C       real    g(2)           : coupling constants                  gvf
C       real    vmass          : mass  of output vector v               
C       real    vwidth         : width of output vector v               
C                                                                       
C output:                                                               
C       complex jio(6)         : vector current          j^mu(<fo|v|fi>)
C
#if defined(CERNLIB_IMPNONE)
      IMPLICIT NONE
#endif
      COMPLEX*16 FI(6),FO(6),JIO(6),C0,C1,C2,C3,CS,D
      REAL*8    G(2),Q(0:3),VMASS,VWIDTH,Q2,VM2,DD
C
      REAL*8 RXZERO, RXONE
      PARAMETER( RXZERO=0.0D0, RXONE=1.0D0 )
      COMPLEX*16 CXIMAG
      LOGICAL FIRST
      SAVE CXIMAG,FIRST
      DATA FIRST/.TRUE./
C
C          Fix compilation with g77
      IF(FIRST) THEN
        FIRST=.FALSE.
        CXIMAG=DCMPLX( RXZERO, RXONE )
      ENDIF
C
      JIO(5) = FO(5)-FI(5)
      JIO(6) = FO(6)-FI(6)
C
      Q(0)=DBLE( JIO(5))
      Q(1)=DBLE( JIO(6))
      Q(2)=DIMAG(JIO(6))
      Q(3)=DIMAG(JIO(5))
      Q2=Q(0)**2-(Q(1)**2+Q(2)**2+Q(3)**2)
      VM2=VMASS**2
C
      IF (VMASS.NE.RXZERO) THEN
C
         D=RXONE/DCMPLX( Q2-VM2 , MAX(SIGN( VMASS*VWIDTH ,Q2),RXZERO) )
C  for the running width, use below instead of the above d.
C      d=r_one/dcmplx( q2-vm2 , max( vwidth*q2/vmass ,r_zero) )
C
         IF (G(2).NE.RXZERO) THEN
C
            C0=  G(1)*( FO(3)*FI(1)+FO(4)*FI(2))
     &          +G(2)*( FO(1)*FI(3)+FO(2)*FI(4))
            C1= -G(1)*( FO(3)*FI(2)+FO(4)*FI(1))
     &          +G(2)*( FO(1)*FI(4)+FO(2)*FI(3))
            C2=( G(1)*( FO(3)*FI(2)-FO(4)*FI(1)) 
     &          +G(2)*(-FO(1)*FI(4)+FO(2)*FI(3)))*CXIMAG
            C3=  G(1)*(-FO(3)*FI(1)+FO(4)*FI(2))
     &          +G(2)*( FO(1)*FI(3)-FO(2)*FI(4))
         ELSE
C
            D=D*G(1)
            C0=  FO(3)*FI(1)+FO(4)*FI(2)
            C1= -FO(3)*FI(2)-FO(4)*FI(1)
            C2=( FO(3)*FI(2)-FO(4)*FI(1))*CXIMAG
            C3= -FO(3)*FI(1)+FO(4)*FI(2)
         END IF
C
         CS=(Q(0)*C0-Q(1)*C1-Q(2)*C2-Q(3)*C3)/VM2
C
         JIO(1) = (C0-CS*Q(0))*D
         JIO(2) = (C1-CS*Q(1))*D
         JIO(3) = (C2-CS*Q(2))*D
         JIO(4) = (C3-CS*Q(3))*D
C
      ELSE
         DD=RXONE/Q2
C
         IF (G(2).NE.RXZERO) THEN
            JIO(1) = ( G(1)*( FO(3)*FI(1)+FO(4)*FI(2))
     &                +G(2)*( FO(1)*FI(3)+FO(2)*FI(4)) )*DD
            JIO(2) = (-G(1)*( FO(3)*FI(2)+FO(4)*FI(1))
     &                +G(2)*( FO(1)*FI(4)+FO(2)*FI(3)) )*DD
            JIO(3) = ( G(1)*( FO(3)*FI(2)-FO(4)*FI(1))
     &                +G(2)*(-FO(1)*FI(4)+FO(2)*FI(3)))
     $           *DCMPLX(RXZERO,DD)
            JIO(4) = ( G(1)*(-FO(3)*FI(1)+FO(4)*FI(2))
     &                +G(2)*( FO(1)*FI(3)-FO(2)*FI(4)) )*DD
C
         ELSE
            DD=DD*G(1)
C
            JIO(1) =  ( FO(3)*FI(1)+FO(4)*FI(2))*DD
            JIO(2) = -( FO(3)*FI(2)+FO(4)*FI(1))*DD
            JIO(3) =  ( FO(3)*FI(2)-FO(4)*FI(1))*DCMPLX(RXZERO,DD)
            JIO(4) =  (-FO(3)*FI(1)+FO(4)*FI(2))*DD
         END IF
      END IF
C
      RETURN
      END
C ----------------------------------------------------------------------
C
      SUBROUTINE JSSXXX(S1,S2,G,VMASS,VWIDTH , JSS)
C
C This subroutine computes an off-shell vector current from the vector- 
C scalar-scalar coupling.  The coupling is absent in the minimal SM in  
C unitary gauge.  The propagator is given in Feynman gauge for a        
C massless vector and in unitary gauge for a massive vector.            
C                                                                       
C INPUT:                                                                
C       complex S1(3)          : first  scalar                        S1
C       complex S2(3)          : second scalar                        S2
C       real    G              : coupling constant (S1 charge)          
C       real    VMASS          : mass  of OUTPUT vector V               
C       real    VWIDTH         : width of OUTPUT vector V               
C                                                                       
C Examples of the coupling constant G for SUSY particles are as follows:
C   -----------------------------------------------------------         
C   |    S1    | (Q,I3) of S1  ||   V=A   |   V=Z   |   V=W   |         
C   -----------------------------------------------------------         
C   | nu~_L    | (  0  , +1/2) ||   ---   |  GZN(1) |  GWF(1) |         
C   | e~_L     | ( -1  , -1/2) ||  GAL(1) |  GZL(1) |  GWF(1) |         
C   | u~_L     | (+2/3 , +1/2) ||  GAU(1) |  GZU(1) |  GWF(1) |         
C   | d~_L     | (-1/3 , -1/2) ||  GAD(1) |  GZD(1) |  GWF(1) |         
C   -----------------------------------------------------------         
C   | e~_R-bar | ( +1  ,  0  ) || -GAL(2) | -GZL(2) | -GWF(2) |         
C   | u~_R-bar | (-2/3 ,  0  ) || -GAU(2) | -GZU(2) | -GWF(2) |         
C   | d~_R-bar | (+1/3 ,  0  ) || -GAD(2) | -GZD(2) | -GWF(2) |         
C   -----------------------------------------------------------         
C where the S1 charge is defined by the flowing-OUT quantum number.     
C                                                                       
C OUTPUT:                                                               
C       complex JSS(6)         : vector current            J^mu(V:S1,S2)
C
#if defined(CERNLIB_IMPNONE)
      IMPLICIT NONE
#endif
      COMPLEX*16 S1(3),S2(3),JSS(6),DG,ADG
      REAL*8  PP(0:3),PA(0:3),Q(0:3),G,VMASS,VWIDTH,Q2,VM2,MP2,MA2,M2D
C
      JSS(5) = S1(2)+S2(2)
      JSS(6) = S1(3)+S2(3)
C
      Q(0)=DBLE( JSS(5))
      Q(1)=DBLE( JSS(6))
      Q(2)=DIMAG(JSS(6))
      Q(3)=DIMAG(JSS(5))
      Q2=Q(0)**2-(Q(1)**2+Q(2)**2+Q(3)**2)
      VM2=VMASS**2
C
      IF (VMASS.EQ.0.) GOTO 10
C
      DG=G/DCMPLX( Q2-VM2, MAX(SIGN( VMASS*VWIDTH ,Q2),0.D0))
C  For the running width, use below instead of the above DG.
C      DG=G/dCMPLX( Q2-VM2 , MAX( VWIDTH*Q2/VMASS ,0.) )
C
      ADG=DG*S1(1)*S2(1)
C
      PP(0)=DBLE( S1(2))
      PP(1)=DBLE( S1(3))
      PP(2)=DIMAG(S1(3))
      PP(3)=DIMAG(S1(2))
      PA(0)=DBLE( S2(2))
      PA(1)=DBLE( S2(3))
      PA(2)=DIMAG(S2(3))
      PA(3)=DIMAG(S2(2))
      MP2=PP(0)**2-(PP(1)**2+PP(2)**2+PP(3)**2)
      MA2=PA(0)**2-(PA(1)**2+PA(2)**2+PA(3)**2)
      M2D=MP2-MA2
C
      JSS(1) = ADG*( (PP(0)-PA(0)) - Q(0)*M2D/VM2)
      JSS(2) = ADG*( (PP(1)-PA(1)) - Q(1)*M2D/VM2)
      JSS(3) = ADG*( (PP(2)-PA(2)) - Q(2)*M2D/VM2)
      JSS(4) = ADG*( (PP(3)-PA(3)) - Q(3)*M2D/VM2)
C
      RETURN
C
  10  ADG=G*S1(1)*S2(1)/Q2
C
      JSS(1) = ADG*DBLE( S1(2)-S2(2))
      JSS(2) = ADG*DBLE( S1(3)-S2(3))
      JSS(3) = ADG*DIMAG(S1(3)-S2(3))
      JSS(4) = ADG*DIMAG(S1(2)-S2(2))
C
      RETURN
      END
C
C
C ----------------------------------------------------------------------
C
      SUBROUTINE JTIOXX(FI,FO,G , JIO)
C
C this subroutine computes an off-shell vector current from an external 
C fermion pair.  the vector boson propagator is not included in this
C routine.
C                                                                       
C input:                                                                
C       complex fi(6)          : flow-in  fermion                   |fi>
C       complex fo(6)          : flow-out fermion                   <fo|
C       real    g(2)           : coupling constants                  gvf
C                                                                       
C output:                                                               
C       complex jio(6)         : vector current          j^mu(<fo|v|fi>)
C
#if defined(CERNLIB_IMPNONE)
      IMPLICIT NONE
#endif
      COMPLEX*16 FI(6),FO(6),JIO(6)
      REAL*8    G(2)
C
      REAL*8 RXZERO, RXONE
      PARAMETER( RXZERO=0.0D0, RXONE=1.0D0 )
      COMPLEX*16 CXIMAG
      LOGICAL FIRST
      SAVE CXIMAG,FIRST
      DATA FIRST/.TRUE./
C
C          Fix compilation with g77
      IF(FIRST) THEN
        FIRST=.FALSE.
        CXIMAG=DCMPLX( RXZERO, RXONE )
      ENDIF
C
      JIO(5) = FO(5)-FI(5)
      JIO(6) = FO(6)-FI(6)
C
      IF ( G(2) .NE. RXZERO ) THEN
         JIO(1) = ( G(1)*( FO(3)*FI(1)+FO(4)*FI(2))
     &             +G(2)*( FO(1)*FI(3)+FO(2)*FI(4)) )
         JIO(2) = (-G(1)*( FO(3)*FI(2)+FO(4)*FI(1))
     &             +G(2)*( FO(1)*FI(4)+FO(2)*FI(3)) )
         JIO(3) = ( G(1)*( FO(3)*FI(2)-FO(4)*FI(1))
     &             +G(2)*(-FO(1)*FI(4)+FO(2)*FI(3)) )*CXIMAG
         JIO(4) = ( G(1)*(-FO(3)*FI(1)+FO(4)*FI(2))
     &             +G(2)*( FO(1)*FI(3)-FO(2)*FI(4)) )
C
      ELSE
         JIO(1) =  ( FO(3)*FI(1)+FO(4)*FI(2))*G(1)
         JIO(2) = -( FO(3)*FI(2)+FO(4)*FI(1))*G(1)
         JIO(3) =  ( FO(3)*FI(2)-FO(4)*FI(1))*DCMPLX(RXZERO,G(1))
         JIO(4) =  (-FO(3)*FI(1)+FO(4)*FI(2))*G(1)
      END IF
C
      RETURN
      END
C ----------------------------------------------------------------------
C
      SUBROUTINE JVSSXX(VC,S1,S2,G,VMASS,VWIDTH , JVSS)
C
C This subroutine computes an off-shell vector current from the vector- 
C vector-scalar-scalar coupling.  The vector propagator is given in     
C Feynman gauge for a massless vector and in unitary gauge for a massive
C vector.                                                               
C                                                                       
C INPUT:                                                                
C       complex VC(6)          : input  vector                        V 
C       complex S1(3)          : first  scalar                        S1
C       complex S2(3)          : second scalar                        S2
C       real    G              : coupling constant                 GVVHH
C       real    VMASS          : mass  of OUTPUT vector V'              
C       real    VWIDTH         : width of OUTPUT vector V'              
C                                                                       
C OUTPUT:                                                               
C       complex JVSS(6)        : vector current         J^mu(V':V,S1,S2)
C
#if defined(CERNLIB_IMPNONE)
      IMPLICIT NONE
#endif
      COMPLEX*16 VC(6),S1(3),S2(3),JVSS(6),DG
      REAL*8    Q(0:3),G,VMASS,VWIDTH,Q2,VK,VM2
C
      JVSS(5) = VC(5)+S1(2)+S2(2)
      JVSS(6) = VC(6)+S1(3)+S2(3)
C
      Q(0)=DBLE( JVSS(5))
      Q(1)=DBLE( JVSS(6))
      Q(2)=DIMAG(JVSS(6))
      Q(3)=DIMAG(JVSS(5))
      Q2=Q(0)**2-(Q(1)**2+Q(2)**2+Q(3)**2)
      VM2=VMASS**2
C
      IF (VMASS.EQ.0.) GOTO 10
C
      DG=G*S1(1)*S2(1)/DCMPLX( Q2-VM2,MAX(SIGN( VMASS*VWIDTH,Q2),0.D0))
C  For the running width, use below instead of the above DG.
C      DG=G*S1(1)*S2(1)/CMPLX( Q2-VM2 , MAX( VWIDTH*Q2/VMASS ,0.))
C
      VK=(Q(0)*VC(1)-Q(1)*VC(2)-Q(2)*VC(3)-Q(3)*VC(4))/VM2
C
      JVSS(1) = DG*(VC(1)-VK*Q(0))
      JVSS(2) = DG*(VC(2)-VK*Q(1))
      JVSS(3) = DG*(VC(3)-VK*Q(2))
      JVSS(4) = DG*(VC(4)-VK*Q(3))
C
      RETURN
C
  10  DG= G*S1(1)*S2(1)/Q2
C
      JVSS(1) = DG*VC(1)
      JVSS(2) = DG*VC(2)
      JVSS(3) = DG*VC(3)
      JVSS(4) = DG*VC(4)
C
      RETURN
      END
C
C
C ----------------------------------------------------------------------
C
      SUBROUTINE JVSXXX(VC,SC,G,VMASS,VWIDTH , JVS)
C
C this subroutine computes an off-shell vector current from the vector- 
C vector-scalar coupling.  the vector propagator is given in feynman    
C gauge for a massless vector and in unitary gauge for a massive vector.
C                                                                       
C input:                                                                
C       complex vc(6)          : input vector                          v
C       complex sc(3)          : input scalar                          s
C       real    g              : coupling constant                  gvvh
C       real    vmass          : mass  of output vector v'              
C       real    vwidth         : width of output vector v'              
C                                                                       
C output:                                                               
C       complex jvs(6)         : vector current             j^mu(v':v,s)
C
#if defined(CERNLIB_IMPNONE)
      IMPLICIT NONE
#endif
      COMPLEX*16 VC(6),SC(3),JVS(6),DG,VK
      REAL*8    Q(0:3),VMASS,VWIDTH,Q2,VM2,G
C
      JVS(5) = VC(5)+SC(2)
      JVS(6) = VC(6)+SC(3)
C
      Q(0)=DBLE( JVS(5))
      Q(1)=DBLE( JVS(6))
      Q(2)=DIMAG(JVS(6))
      Q(3)=DIMAG(JVS(5))
      Q2=Q(0)**2-(Q(1)**2+Q(2)**2+Q(3)**2)
      VM2=VMASS**2
C
      IF (VMASS.EQ.0.) GOTO 10
C
      DG=G*SC(1)/DCMPLX( Q2-VM2 , MAX(DSIGN( VMASS*VWIDTH ,Q2),0.D0) )
C  for the running width, use below instead of the above dg.
C      dg=g*sc(1)/dcmplx( q2-vm2 , max( vwidth*q2/vmass ,0.) )
C
      VK=(-Q(0)*VC(1)+Q(1)*VC(2)+Q(2)*VC(3)+Q(3)*VC(4))/VM2
C
      JVS(1) = DG*(Q(0)*VK+VC(1))
      JVS(2) = DG*(Q(1)*VK+VC(2))
      JVS(3) = DG*(Q(2)*VK+VC(3))
      JVS(4) = DG*(Q(3)*VK+VC(4))
C
      RETURN
C
  10  DG=G*SC(1)/Q2
C
      JVS(1) = DG*VC(1)
      JVS(2) = DG*VC(2)
      JVS(3) = DG*VC(3)
      JVS(4) = DG*VC(4)
C
      RETURN
      END


C
C ----------------------------------------------------------------------
C
      SUBROUTINE JVVXXX(V1,V2,G,VMASS,VWIDTH , JVV)
C
C this subroutine computes an off-shell vector current from the three-  
C point gauge boson coupling.  the vector propagator is given in feynman
C gauge for a massless vector and in unitary gauge for a massive vector.
C                                                                       
C input:                                                                
C       complex v1(6)          : first  vector                        v1
C       complex v2(6)          : second vector                        v2
C       real    g              : coupling constant (see the table below)
C       real    vmass          : mass  of output vector v               
C       real    vwidth         : width of output vector v               
C                                                                       
C the possible sets of the inputs are as follows:                       
C    ------------------------------------------------------------------ 
C    |   v1   |   v2   |  jvv   |      g       |   vmass  |  vwidth   | 
C    ------------------------------------------------------------------ 
C    |   w-   |   w+   |  a/z   |  gwwa/gwwz   | 0./zmass | 0./zwidth | 
C    | w3/a/z |   w-   |  w+    | gw/gwwa/gwwz |   wmass  |  wwidth   | 
C    |   w+   | w3/a/z |  w-    | gw/gwwa/gwwz |   wmass  |  wwidth   | 
C    ------------------------------------------------------------------ 
C where all the bosons are defined by the flowing-out quantum number.   
C                                                                       
C output:                                                               
C       complex jvv(6)         : vector current            j^mu(v:v1,v2)
C
#if defined(CERNLIB_IMPNONE)
      IMPLICIT NONE
#endif
      COMPLEX*16 V1(6),V2(6),JVV(6),J12(0:3),JS,DG,
     &        SV1,SV2,S11,S12,S21,S22,V12
      REAL*8    P1(0:3),P2(0:3),Q(0:3),G,VMASS,VWIDTH,GS,S,VM2,M1,M2
C
      REAL*8 RXZERO
      PARAMETER( RXZERO=0.0D0 )
C
      JVV(5) = V1(5)+V2(5)
      JVV(6) = V1(6)+V2(6)
C
      P1(0)=DBLE( V1(5))
      P1(1)=DBLE( V1(6))
      P1(2)=DIMAG(V1(6))
      P1(3)=DIMAG(V1(5))
      P2(0)=DBLE( V2(5))
      P2(1)=DBLE( V2(6))
      P2(2)=DIMAG(V2(6))
      P2(3)=DIMAG(V2(5))
      Q(0)=-DBLE( JVV(5))
      Q(1)=-DBLE( JVV(6))
      Q(2)=-DIMAG(JVV(6))
      Q(3)=-DIMAG(JVV(5))
      S=Q(0)**2-(Q(1)**2+Q(2)**2+Q(3)**2)
C
      V12=V1(1)*V2(1)-V1(2)*V2(2)-V1(3)*V2(3)-V1(4)*V2(4)
      SV1= (P2(0)-Q(0))*V1(1) -(P2(1)-Q(1))*V1(2)
     &    -(P2(2)-Q(2))*V1(3) -(P2(3)-Q(3))*V1(4)
      SV2=-(P1(0)-Q(0))*V2(1) +(P1(1)-Q(1))*V2(2)
     &    +(P1(2)-Q(2))*V2(3) +(P1(3)-Q(3))*V2(4)
      J12(0)=(P1(0)-P2(0))*V12 +SV1*V2(1) +SV2*V1(1)
      J12(1)=(P1(1)-P2(1))*V12 +SV1*V2(2) +SV2*V1(2)
      J12(2)=(P1(2)-P2(2))*V12 +SV1*V2(3) +SV2*V1(3)
      J12(3)=(P1(3)-P2(3))*V12 +SV1*V2(4) +SV2*V1(4)
C
      IF ( VMASS .NE. RXZERO ) THEN
         VM2=VMASS**2
         M1=P1(0)**2-(P1(1)**2+P1(2)**2+P1(3)**2)
         M2=P2(0)**2-(P2(1)**2+P2(2)**2+P2(3)**2)
         S11=P1(0)*V1(1)-P1(1)*V1(2)-P1(2)*V1(3)-P1(3)*V1(4)
         S12=P1(0)*V2(1)-P1(1)*V2(2)-P1(2)*V2(3)-P1(3)*V2(4)
         S21=P2(0)*V1(1)-P2(1)*V1(2)-P2(2)*V1(3)-P2(3)*V1(4)
         S22=P2(0)*V2(1)-P2(1)*V2(2)-P2(2)*V2(3)-P2(3)*V2(4)
         JS=(V12*(-M1+M2) +S11*S12 -S21*S22)/VM2
C
         DG=-G/DCMPLX( S-VM2 , MAX(SIGN( VMASS*VWIDTH ,S),RXZERO) )
C
C  for the running width, use below instead of the above dg.
C         dg=-g/dcmplx( s-vm2 , max( vwidth*s/vmass ,r_zero) )
C
         JVV(1) = DG*(J12(0)-Q(0)*JS)
         JVV(2) = DG*(J12(1)-Q(1)*JS)
         JVV(3) = DG*(J12(2)-Q(2)*JS)
         JVV(4) = DG*(J12(3)-Q(3)*JS)
C
      ELSE
         GS=-G/S
C
         JVV(1) = GS*J12(0)
         JVV(2) = GS*J12(1)
         JVV(3) = GS*J12(2)
         JVV(4) = GS*J12(3)
      END IF
C
      RETURN
      END
C
C ----------------------------------------------------------------------
C
      SUBROUTINE JW3WXX(W1,W2,W3,G1,G2,WMASS,WWIDTH,VMASS,VWIDTH , JW3W)
C
C this subroutine computes an off-shell w+, w-, w3, z or photon current 
C from the four-point gauge boson coupling, including the contributions 
C of w exchange diagrams.  the vector propagator is given in feynman    
C gauge for a photon and in unitary gauge for w and z bosons.  if one   
C sets wmass=0.0, then the ggg-->g current is given (see sect 2.9.1 of 
C the manual).                                                          
C                                                                       
C input:                                                                
C       complex w1(6)          : first  vector                        w1
C       complex w2(6)          : second vector                        w2
C       complex w3(6)          : third  vector                        w3
C       real    g1             : first  coupling constant               
C       real    g2             : second coupling constant               
C                                                  (see the table below)
C       real    wmass          : mass  of internal w                    
C       real    wwidth         : width of internal w                    
C       real    vmass          : mass  of output w'                     
C       real    vwidth         : width of output w'                     
C                                                                       
C the possible sets of the inputs are as follows:                       
C   ------------------------------------------------------------------- 
C   |  w1  |  w2  |  w3  | g1 | g2 |wmass|wwidth|vmass|vwidth || jw3w | 
C   ------------------------------------------------------------------- 
C   |  w-  |  w3  |  w+  | gw |gwwz|wmass|wwidth|zmass|zwidth ||  z   | 
C   |  w-  |  w3  |  w+  | gw |gwwa|wmass|wwidth|  0. |  0.   ||  a   | 
C   |  w-  |  z   |  w+  |gwwz|gwwz|wmass|wwidth|zmass|zwidth ||  z   | 
C   |  w-  |  z   |  w+  |gwwz|gwwa|wmass|wwidth|  0. |  0.   ||  a   | 
C   |  w-  |  a   |  w+  |gwwa|gwwz|wmass|wwidth|zmass|zwidth ||  z   | 
C   |  w-  |  a   |  w+  |gwwa|gwwa|wmass|wwidth|  0. |  0.   ||  a   | 
C   ------------------------------------------------------------------- 
C   |  w3  |  w-  |  w3  | gw | gw |wmass|wwidth|wmass|wwidth ||  w+  | 
C   |  w3  |  w+  |  w3  | gw | gw |wmass|wwidth|wmass|wwidth ||  w-  | 
C   |  w3  |  w-  |  z   | gw |gwwz|wmass|wwidth|wmass|wwidth ||  w+  | 
C   |  w3  |  w+  |  z   | gw |gwwz|wmass|wwidth|wmass|wwidth ||  w-  | 
C   |  w3  |  w-  |  a   | gw |gwwa|wmass|wwidth|wmass|wwidth ||  w+  | 
C   |  w3  |  w+  |  a   | gw |gwwa|wmass|wwidth|wmass|wwidth ||  w-  | 
C   |  z   |  w-  |  z   |gwwz|gwwz|wmass|wwidth|wmass|wwidth ||  w+  | 
C   |  z   |  w+  |  z   |gwwz|gwwz|wmass|wwidth|wmass|wwidth ||  w-  | 
C   |  z   |  w-  |  a   |gwwz|gwwa|wmass|wwidth|wmass|wwidth ||  w+  | 
C   |  z   |  w+  |  a   |gwwz|gwwa|wmass|wwidth|wmass|wwidth ||  w-  | 
C   |  a   |  w-  |  a   |gwwa|gwwa|wmass|wwidth|wmass|wwidth ||  w+  | 
C   |  a   |  w+  |  a   |gwwa|gwwa|wmass|wwidth|wmass|wwidth ||  w-  | 
C   ------------------------------------------------------------------- 
C where all the bosons are defined by the flowing-out quantum number.   
C                                                                       
C output:                                                               
C       complex jw3w(6)        : w current             j^mu(w':w1,w2,w3)
C
#if defined(CERNLIB_IMPNONE)
      IMPLICIT NONE
#endif
      COMPLEX*16  W1(6),W2(6),W3(6),JW3W(6)
      COMPLEX*16 DW1(0:3),DW2(0:3),DW3(0:3),
     &           JJ(0:3),J4(0:3),
     &           DV,W12,W32,W13,
     &           JQ
      REAL*8     G1,G2,WMASS,WWIDTH,VMASS,VWIDTH
      REAL*8     P1(0:3),P2(0:3),P3(0:3),Q(0:3),
     &           DG2,DMV,DWV,MV2,Q2
C
      REAL*8 RXZERO
      PARAMETER( RXZERO=0.0D0 )
C
      JW3W(5) = W1(5)+W2(5)+W3(5)
      JW3W(6) = W1(6)+W2(6)+W3(6)
C
      DW1(0)=DCMPLX(W1(1))
      DW1(1)=DCMPLX(W1(2))
      DW1(2)=DCMPLX(W1(3))
      DW1(3)=DCMPLX(W1(4))
      DW2(0)=DCMPLX(W2(1))
      DW2(1)=DCMPLX(W2(2))
      DW2(2)=DCMPLX(W2(3))
      DW2(3)=DCMPLX(W2(4))
      DW3(0)=DCMPLX(W3(1))
      DW3(1)=DCMPLX(W3(2))
      DW3(2)=DCMPLX(W3(3))
      DW3(3)=DCMPLX(W3(4))
      P1(0)=DBLE(      W1(5))
      P1(1)=DBLE(      W1(6))
      P1(2)=DBLE(DIMAG(W1(6)))
      P1(3)=DBLE(DIMAG(W1(5)))
      P2(0)=DBLE(      W2(5))
      P2(1)=DBLE(      W2(6))
      P2(2)=DBLE(DIMAG(W2(6)))
      P2(3)=DBLE(DIMAG(W2(5)))
      P3(0)=DBLE(      W3(5))
      P3(1)=DBLE(      W3(6))
      P3(2)=DBLE(DIMAG(W3(6)))
      P3(3)=DBLE(DIMAG(W3(5)))
      Q(0)=-(P1(0)+P2(0)+P3(0))
      Q(1)=-(P1(1)+P2(1)+P3(1))
      Q(2)=-(P1(2)+P2(2)+P3(2))
      Q(3)=-(P1(3)+P2(3)+P3(3))


      Q2 =Q(0)**2 -(Q(1)**2 +Q(2)**2 +Q(3)**2)
      DG2=DBLE(G1)*DBLE(G2)
      DMV=DBLE(VMASS)
      DWV=DBLE(VWIDTH)
      MV2=DMV**2
      IF (VMASS.EQ. RXZERO) THEN
      DV = 1.0D0/DCMPLX( Q2 )
      ELSE
      DV = 1.0D0/DCMPLX( Q2 -MV2 , DMAX1(DSIGN(DMV*DWV,Q2 ),0.D0) )
      ENDIF
C  for the running width, use below instead of the above dv.
C      dv = 1.0d0/dcmplx( q2 -mv2 , dmax1(dwv*q2/dmv,0.d0) )
C
      W12=DW1(0)*DW2(0)-DW1(1)*DW2(1)-DW1(2)*DW2(2)-DW1(3)*DW2(3)
      W32=DW3(0)*DW2(0)-DW3(1)*DW2(1)-DW3(2)*DW2(2)-DW3(3)*DW2(3)
C
      IF ( WMASS .NE. RXZERO ) THEN
         W13=DW1(0)*DW3(0)-DW1(1)*DW3(1)-DW1(2)*DW3(2)-DW1(3)*DW3(3)
C
         J4(0)=DG2*( DW1(0)*W32 + DW3(0)*W12 - 2.D0*DW2(0)*W13 )
         J4(1)=DG2*( DW1(1)*W32 + DW3(1)*W12 - 2.D0*DW2(1)*W13 )
         J4(2)=DG2*( DW1(2)*W32 + DW3(2)*W12 - 2.D0*DW2(2)*W13 )
         J4(3)=DG2*( DW1(3)*W32 + DW3(3)*W12 - 2.D0*DW2(3)*W13 )
C
         JJ(0)=J4(0)
         JJ(1)=J4(1)
         JJ(2)=J4(2)
         JJ(3)=J4(3)

      ELSE
C
         W12=DW1(0)*DW2(0)-DW1(1)*DW2(1)-DW1(2)*DW2(2)-DW1(3)*DW2(3)
         W32=DW3(0)*DW2(0)-DW3(1)*DW2(1)-DW3(2)*DW2(2)-DW3(3)*DW2(3)
         W13=DW1(0)*DW3(0)-DW1(1)*DW3(1)-DW1(2)*DW3(2)-DW1(3)*DW3(3)
C
         J4(0)=DG2*( DW1(0)*W32 - DW2(0)*W13 )
         J4(1)=DG2*( DW1(1)*W32 - DW2(1)*W13 )
         J4(2)=DG2*( DW1(2)*W32 - DW2(2)*W13 )
         J4(3)=DG2*( DW1(3)*W32 - DW2(3)*W13 )
C
         JJ(0)=J4(0)
         JJ(1)=J4(1)
         JJ(2)=J4(2)
         JJ(3)=J4(3)

      END IF
C
      IF ( VMASS .NE. RXZERO ) THEN
C
         JQ=(JJ(0)*Q(0)-JJ(1)*Q(1)-JJ(2)*Q(2)-JJ(3)*Q(3))/MV2
C
         JW3W(1) = DCMPLX( (JJ(0)-JQ*Q(0))*DV )
         JW3W(2) = DCMPLX( (JJ(1)-JQ*Q(1))*DV )
         JW3W(3) = DCMPLX( (JJ(2)-JQ*Q(2))*DV )
         JW3W(4) = DCMPLX( (JJ(3)-JQ*Q(3))*DV )
C
      ELSE
C
         JW3W(1) = DCMPLX( JJ(0)*DV )
         JW3W(2) = DCMPLX( JJ(1)*DV )
         JW3W(3) = DCMPLX( JJ(2)*DV )
         JW3W(4) = DCMPLX( JJ(3)*DV )
      END IF
C
      RETURN
      END
C
C ----------------------------------------------------------------------
C
      SUBROUTINE JWWWXX(W1,W2,W3,GWWA,GWWZ,ZMASS,ZWIDTH,WMASS,WWIDTH ,
     &                  JWWW)
C
C this subroutine computes an off-shell w+/w- current from the four-    
C point gauge boson coupling, including the contributions of photon and 
C z exchanges.  the vector propagators for the output w and the internal
C z bosons are given in unitary gauge, and that of the internal photon  
C is given in feynman gauge.                                            
C                                                                       
C input:                                                                
C       complex w1(6)          : first  vector                        w1
C       complex w2(6)          : second vector                        w2
C       complex w3(6)          : third  vector                        w3
C       real    gwwa           : coupling constant of w and a       gwwa
C       real    gwwz           : coupling constant of w and z       gwwz
C       real    zmass          : mass  of internal z                    
C       real    zwidth         : width of internal z                    
C       real    wmass          : mass  of output w                      
C       real    wwidth         : width of output w                      
C                                                                       
C the possible sets of the inputs are as follows:                       
C   ------------------------------------------------------------------- 
C   |  w1  |  w2  |  w3  |gwwa|gwwz|zmass|zwidth|wmass|wwidth || jwww | 
C   ------------------------------------------------------------------- 
C   |  w-  |  w+  |  w-  |gwwa|gwwz|zmass|zwidth|wmass|wwidth ||  w+  | 
C   |  w+  |  w-  |  w+  |gwwa|gwwz|zmass|zwidth|wmass|wwidth ||  w-  | 
C   ------------------------------------------------------------------- 
C where all the bosons are defined by the flowing-out quantum number.   
C                                                                       
C output:                                                               
C       complex jwww(6)        : w current             j^mu(w':w1,w2,w3)
C
#if defined(CERNLIB_IMPNONE)
      IMPLICIT NONE
#endif
      COMPLEX*16  W1(6),W2(6),W3(6),JWWW(6)
      COMPLEX*16 DW1(0:3),DW2(0:3),DW3(0:3),
     &           JJ(0:3),JS(0:3),JT(0:3),J4(0:3),
     &           JT12(0:3),JT32(0:3),J12(0:3),J32(0:3),
     &           DZS,DZT,DW,W12,W32,W13,P1W2,P2W1,P3W2,P2W3,
     &           JK12,JK32,JSW3,JTW1,P3JS,KSW3,P1JT,KTW1,JQ
      REAL*8     GWWA,GWWZ,ZMASS,ZWIDTH,WMASS,WWIDTH
      REAL*8     P1(0:3),P2(0:3),P3(0:3),Q(0:3),KS(0:3),KT(0:3),
     &           DGWWA2,DGWWZ2,DGW2,DMZ,DWZ,DMW,DWW,MZ2,MW2,Q2,KS2,KT2,
     &           DAS,DAT
C
      JWWW(5) = W1(5)+W2(5)+W3(5)
      JWWW(6) = W1(6)+W2(6)+W3(6)
C
      DW1(0)=DCMPLX(W1(1))
      DW1(1)=DCMPLX(W1(2))
      DW1(2)=DCMPLX(W1(3))
      DW1(3)=DCMPLX(W1(4))
      DW2(0)=DCMPLX(W2(1))
      DW2(1)=DCMPLX(W2(2))
      DW2(2)=DCMPLX(W2(3))
      DW2(3)=DCMPLX(W2(4))
      DW3(0)=DCMPLX(W3(1))
      DW3(1)=DCMPLX(W3(2))
      DW3(2)=DCMPLX(W3(3))
      DW3(3)=DCMPLX(W3(4))
      P1(0)=DBLE(      W1(5))
      P1(1)=DBLE(      W1(6))
      P1(2)=DBLE(DIMAG(W1(6)))
      P1(3)=DBLE(DIMAG(W1(5)))
      P2(0)=DBLE(      W2(5))
      P2(1)=DBLE(      W2(6))
      P2(2)=DBLE(DIMAG(W2(6)))
      P2(3)=DBLE(DIMAG(W2(5)))
      P3(0)=DBLE(      W3(5))
      P3(1)=DBLE(      W3(6))
      P3(2)=DBLE(DIMAG(W3(6)))
      P3(3)=DBLE(DIMAG(W3(5)))
      Q(0)=-(P1(0)+P2(0)+P3(0))
      Q(1)=-(P1(1)+P2(1)+P3(1))
      Q(2)=-(P1(2)+P2(2)+P3(2))
      Q(3)=-(P1(3)+P2(3)+P3(3))
      KS(0)=P1(0)+P2(0)
      KS(1)=P1(1)+P2(1)
      KS(2)=P1(2)+P2(2)
      KS(3)=P1(3)+P2(3)
      KT(0)=P2(0)+P3(0)
      KT(1)=P2(1)+P3(1)
      KT(2)=P2(2)+P3(2)
      KT(3)=P2(3)+P3(3)
      Q2 =Q(0)**2 -(Q(1)**2 +Q(2)**2 +Q(3)**2)
      KS2=KS(0)**2-(KS(1)**2+KS(2)**2+KS(3)**2)
      KT2=KT(0)**2-(KT(1)**2+KT(2)**2+KT(3)**2)
      DGWWA2=DBLE(GWWA)**2
      DGWWZ2=DBLE(GWWZ)**2
      DGW2  =DGWWA2+DGWWZ2
      DMZ=DBLE(ZMASS)
      DWZ=DBLE(ZWIDTH)
      DMW=DBLE(WMASS)
      DWW=DBLE(WWIDTH)
      MZ2=DMZ**2
      MW2=DMW**2
C
      DAS=-DGWWA2/KS2
      DAT=-DGWWA2/KT2
      DZS=-DGWWZ2/DCMPLX( KS2-MZ2 , DMAX1(DSIGN(DMZ*DWZ,KS2),0.D0) )
      DZT=-DGWWZ2/DCMPLX( KT2-MZ2 , DMAX1(DSIGN(DMZ*DWZ,KT2),0.D0) )
      DW =-1.0D0/DCMPLX( Q2 -MW2 , DMAX1(DSIGN(DMW*DWW,Q2 ),0.D0) )
C  for the running width, use below instead of the above dw.
C      dw =-1.0d0/dcmplx( q2 -mw2 , dmax1(dww*q2/dmw,0.d0) )
C
      W12=DW1(0)*DW2(0)-DW1(1)*DW2(1)-DW1(2)*DW2(2)-DW1(3)*DW2(3)
      W32=DW3(0)*DW2(0)-DW3(1)*DW2(1)-DW3(2)*DW2(2)-DW3(3)*DW2(3)
C
      P1W2= (P1(0)+KS(0))*DW2(0)-(P1(1)+KS(1))*DW2(1)
     &     -(P1(2)+KS(2))*DW2(2)-(P1(3)+KS(3))*DW2(3)
      P2W1= (P2(0)+KS(0))*DW1(0)-(P2(1)+KS(1))*DW1(1)
     &     -(P2(2)+KS(2))*DW1(2)-(P2(3)+KS(3))*DW1(3)
      P3W2= (P3(0)+KT(0))*DW2(0)-(P3(1)+KT(1))*DW2(1)
     &     -(P3(2)+KT(2))*DW2(2)-(P3(3)+KT(3))*DW2(3)
      P2W3= (P2(0)+KT(0))*DW3(0)-(P2(1)+KT(1))*DW3(1)
     &     -(P2(2)+KT(2))*DW3(2)-(P2(3)+KT(3))*DW3(3)
C
      JT12(0)= (P1(0)-P2(0))*W12 + P2W1*DW2(0) - P1W2*DW1(0)
      JT12(1)= (P1(1)-P2(1))*W12 + P2W1*DW2(1) - P1W2*DW1(1)
      JT12(2)= (P1(2)-P2(2))*W12 + P2W1*DW2(2) - P1W2*DW1(2)
      JT12(3)= (P1(3)-P2(3))*W12 + P2W1*DW2(3) - P1W2*DW1(3)
      JT32(0)= (P3(0)-P2(0))*W32 + P2W3*DW2(0) - P3W2*DW3(0)
      JT32(1)= (P3(1)-P2(1))*W32 + P2W3*DW2(1) - P3W2*DW3(1)
      JT32(2)= (P3(2)-P2(2))*W32 + P2W3*DW2(2) - P3W2*DW3(2)
      JT32(3)= (P3(3)-P2(3))*W32 + P2W3*DW2(3) - P3W2*DW3(3)
C
      JK12=(JT12(0)*KS(0)-JT12(1)*KS(1)-JT12(2)*KS(2)-JT12(3)*KS(3))/MZ2
      JK32=(JT32(0)*KT(0)-JT32(1)*KT(1)-JT32(2)*KT(2)-JT32(3)*KT(3))/MZ2
C
      J12(0)=JT12(0)*(DAS+DZS)-KS(0)*JK12*DZS
      J12(1)=JT12(1)*(DAS+DZS)-KS(1)*JK12*DZS
      J12(2)=JT12(2)*(DAS+DZS)-KS(2)*JK12*DZS
      J12(3)=JT12(3)*(DAS+DZS)-KS(3)*JK12*DZS
      J32(0)=JT32(0)*(DAT+DZT)-KT(0)*JK32*DZT
      J32(1)=JT32(1)*(DAT+DZT)-KT(1)*JK32*DZT
      J32(2)=JT32(2)*(DAT+DZT)-KT(2)*JK32*DZT
      J32(3)=JT32(3)*(DAT+DZT)-KT(3)*JK32*DZT
C
      JSW3=J12(0)*DW3(0)-J12(1)*DW3(1)-J12(2)*DW3(2)-J12(3)*DW3(3)
      JTW1=J32(0)*DW1(0)-J32(1)*DW1(1)-J32(2)*DW1(2)-J32(3)*DW1(3)
C
      P3JS= (P3(0)-Q(0))*J12(0)-(P3(1)-Q(1))*J12(1)
     &     -(P3(2)-Q(2))*J12(2)-(P3(3)-Q(3))*J12(3)
      KSW3= (KS(0)-Q(0))*DW3(0)-(KS(1)-Q(1))*DW3(1)
     &     -(KS(2)-Q(2))*DW3(2)-(KS(3)-Q(3))*DW3(3)
      P1JT= (P1(0)-Q(0))*J32(0)-(P1(1)-Q(1))*J32(1)
     &     -(P1(2)-Q(2))*J32(2)-(P1(3)-Q(3))*J32(3)
      KTW1= (KT(0)-Q(0))*DW1(0)-(KT(1)-Q(1))*DW1(1)
     &     -(KT(2)-Q(2))*DW1(2)-(KT(3)-Q(3))*DW1(3)
C
      JS(0)= (KS(0)-P3(0))*JSW3 + P3JS*DW3(0) - KSW3*J12(0)
      JS(1)= (KS(1)-P3(1))*JSW3 + P3JS*DW3(1) - KSW3*J12(1)
      JS(2)= (KS(2)-P3(2))*JSW3 + P3JS*DW3(2) - KSW3*J12(2)
      JS(3)= (KS(3)-P3(3))*JSW3 + P3JS*DW3(3) - KSW3*J12(3)
      JT(0)= (KT(0)-P1(0))*JTW1 + P1JT*DW1(0) - KTW1*J32(0)
      JT(1)= (KT(1)-P1(1))*JTW1 + P1JT*DW1(1) - KTW1*J32(1)
      JT(2)= (KT(2)-P1(2))*JTW1 + P1JT*DW1(2) - KTW1*J32(2)
      JT(3)= (KT(3)-P1(3))*JTW1 + P1JT*DW1(3) - KTW1*J32(3)
C
      W13=DW1(0)*DW3(0)-DW1(1)*DW3(1)-DW1(2)*DW3(2)-DW1(3)*DW3(3)
C
      J4(0)=DGW2*( DW1(0)*W32 + DW3(0)*W12 - 2.D0*DW2(0)*W13 )
      J4(1)=DGW2*( DW1(1)*W32 + DW3(1)*W12 - 2.D0*DW2(1)*W13 )
      J4(2)=DGW2*( DW1(2)*W32 + DW3(2)*W12 - 2.D0*DW2(2)*W13 )
      J4(3)=DGW2*( DW1(3)*W32 + DW3(3)*W12 - 2.D0*DW2(3)*W13 )
C
C      jj(0)=js(0)+jt(0)+j4(0)
C      jj(1)=js(1)+jt(1)+j4(1)
C      jj(2)=js(2)+jt(2)+j4(2)
C      jj(3)=js(3)+jt(3)+j4(3)

      JJ(0)=J4(0)
      JJ(1)=J4(1)
      JJ(2)=J4(2)
      JJ(3)=J4(3)
C
      JQ=(JJ(0)*Q(0)-JJ(1)*Q(1)-JJ(2)*Q(2)-JJ(3)*Q(3))/MW2
C

      JWWW(1) = DCMPLX( (JJ(0)-JQ*Q(0))*DW )
      JWWW(2) = DCMPLX( (JJ(1)-JQ*Q(1))*DW )
      JWWW(3) = DCMPLX( (JJ(2)-JQ*Q(2))*DW )
      JWWW(4) = DCMPLX( (JJ(3)-JQ*Q(3))*DW )
C
      RETURN
      END

C
C ----------------------------------------------------------------------
C
      SUBROUTINE MOM2CX(ESUM,MASS1,MASS2,COSTH1,PHI1 , P1,P2)
C
C This subroutine sets up two four-momenta in the two particle rest     
C frame.                                                                
C                                                                       
C INPUT:                                                                
C       real    ESUM           : energy sum of particle 1 and 2         
C       real    MASS1          : mass            of particle 1          
C       real    MASS2          : mass            of particle 2          
C       real    COSTH1         : cos(theta)      of particle 1          
C       real    PHI1           : azimuthal angle of particle 1          
C                                                                       
C OUTPUT:                                                               
C       real    P1(0:3)        : four-momentum of particle 1            
C       real    P2(0:3)        : four-momentum of particle 2            
C
#if defined(CERNLIB_IMPNONE)
      IMPLICIT NONE
#endif
      REAL*8    P1(0:3),P2(0:3),
     &        ESUM,MASS1,MASS2,COSTH1,PHI1,MD2,ED,PP,SINTH1
C
      MD2=(MASS1-MASS2)*(MASS1+MASS2)
      ED=MD2/ESUM
      IF (MASS1*MASS2.EQ.0.) THEN
      PP=(ESUM-ABS(ED))*0.5D0
C
      ELSE
      PP=SQRT((MD2/ESUM)**2-2.0D0*(MASS1**2+MASS2**2)+ESUM**2)*0.5D0
      ENDIF
      SINTH1=SQRT((1.0D0-COSTH1)*(1.0D0+COSTH1))
C
      P1(0) = MAX((ESUM+ED)*0.5D0,0.D0)
      P1(1) = PP*SINTH1*COS(PHI1)
      P1(2) = PP*SINTH1*SIN(PHI1)
      P1(3) = PP*COSTH1
C
      P2(0) = MAX((ESUM-ED)*0.5D0,0.D0)
      P2(1) = -P1(1)
      P2(2) = -P1(2)
      P2(3) = -P1(3)
C
      RETURN
      END
C **********************************************************************
C
      SUBROUTINE MOMNTX(ENERGY,MASS,COSTH,PHI , P)
C
C This subroutine sets up a four-momentum from the four inputs.         
C                                                                       
C INPUT:                                                                
C       real    ENERGY         : energy                                 
C       real    MASS           : mass                                   
C       real    COSTH          : cos(theta)                             
C       real    PHI            : azimuthal angle                        
C                                                                       
C OUTPUT:                                                               
C       real    P(0:3)         : four-momentum                          
C
#if defined(CERNLIB_IMPNONE)
      IMPLICIT NONE
#endif
      REAL*8    P(0:3),ENERGY,MASS,COSTH,PHI,PP,SINTH
C
      P(0) = ENERGY
      IF (ENERGY.EQ.MASS) THEN
         P(1) = 0.
         P(2) = 0.
         P(3) = 0.
      ELSE
         PP=SQRT((ENERGY-MASS)*(ENERGY+MASS))
         SINTH=SQRT((1.-COSTH)*(1.+COSTH))
         P(3) = PP*COSTH
         IF (PHI.EQ.0.) THEN
            P(1) = PP*SINTH
            P(2) = 0.
         ELSE
            P(1) = PP*SINTH*COS(PHI)
            P(2) = PP*SINTH*SIN(PHI)
         ENDIF
      ENDIF
      RETURN
      END
C
C
C
C	Subroutine returns the desired fermion or
C	anti-fermion anti-spinor. ie., <f|
C	A replacement for the HELAS routine OXXXXX
C
C	Adam Duff,  1992 August 31
C	<duff@phenom.physics.wisc.edu>
C
      SUBROUTINE OXXXXX(P,FMASS,NHEL,NSF,FO)
C
C          P            IN: FOUR VECTOR MOMENTUM
C          FMASS        IN: FERMION MASS
C          NHEL         IN: ANTI-SPINOR HELICITY, -1 OR 1
C          NSF          IN: -1=ANTIFERMION, 1=FERMION
C          FO           OUT: FERMION WAVEFUNCTION
C
C declare input/output variables
C
#if defined(CERNLIB_IMPNONE)
      IMPLICIT NONE
#endif
      COMPLEX*16 FO(6)
      INTEGER*4 NHEL, NSF
      REAL*8 P(0:3), FMASS
C
C declare local variables
C
      REAL*8 RXZERO, RXONE, RXTWO
      PARAMETER( RXZERO=0.0D0, RXONE=1.0D0, RXTWO=2.0D0 )
      REAL*8 PLAT, PABS, OMEGAP, OMEGAM, RS2PA, SPAZ
      COMPLEX*16 CXZERO
      LOGICAL FIRST
      SAVE CXZERO,FIRST
      DATA FIRST/.TRUE./
C
C          Fix compilation with g77
      IF(FIRST) THEN
        FIRST=.FALSE.
        CXZERO=DCMPLX( RXZERO, RXZERO )
      ENDIF
C
C define kinematic parameters
C
      FO(5) = DCMPLX( P(0), P(3) ) * NSF
      FO(6) = DCMPLX( P(1), P(2) ) * NSF
      PLAT = SQRT( P(1)**2 + P(2)**2 )
      PABS = SQRT( P(1)**2 + P(2)**2 + P(3)**2 )
      OMEGAP = SQRT( P(0) + PABS )
C
C do massive fermion case
C
      IF ( FMASS .NE. RXZERO ) THEN
         OMEGAM = FMASS / OMEGAP
         IF ( NSF .EQ. 1 ) THEN
            IF ( NHEL .EQ. 1 ) THEN
               IF ( P(3) .GE. RXZERO ) THEN
                  IF ( PLAT .EQ. RXZERO ) THEN
                     FO(1) = DCMPLX( OMEGAP, RXZERO )
                     FO(2) = CXZERO
                     FO(3) = DCMPLX( OMEGAM, RXZERO )
                     FO(4) = CXZERO
                  ELSE
                     RS2PA = RXONE / SQRT( RXTWO * PABS )
                     SPAZ = SQRT( PABS + P(3) )
                     FO(1) = OMEGAP * RS2PA
     &                     * DCMPLX( SPAZ, RXZERO )
                     FO(2) = OMEGAP * RS2PA / SPAZ
     &                     * DCMPLX( P(1), -P(2) )
                     FO(3) = OMEGAM * RS2PA
     &                     * DCMPLX( SPAZ, RXZERO )
                     FO(4) = OMEGAM * RS2PA / SPAZ
     &                     * DCMPLX( P(1), -P(2) )
                  END IF
               ELSE
                  IF ( PLAT .EQ. RXZERO ) THEN
                     FO(1) = CXZERO
                     FO(2) = DCMPLX( OMEGAP, RXZERO )
                     FO(3) = CXZERO
                     FO(4) = DCMPLX( OMEGAM, RXZERO )
                  ELSE
                     RS2PA = RXONE / SQRT( RXTWO * PABS )
                     SPAZ = SQRT( PABS - P(3) )
                     FO(1) = OMEGAP * RS2PA / SPAZ
     &                     * DCMPLX( PLAT, RXZERO )
                     FO(2) = OMEGAP * RS2PA * SPAZ / PLAT
     &                     * DCMPLX( P(1), -P(2) )
                     FO(3) = OMEGAM * RS2PA / SPAZ
     &                     * DCMPLX( PLAT, RXZERO )
                     FO(4) = OMEGAM * RS2PA * SPAZ / PLAT
     &                     * DCMPLX( P(1), -P(2) )
                  END IF
               END IF
            ELSE IF ( NHEL .EQ. -1 ) THEN
               IF ( P(3) .GE. RXZERO ) THEN
                  IF ( PLAT .EQ. RXZERO ) THEN
                     FO(1) = CXZERO
                     FO(2) = DCMPLX( OMEGAM, RXZERO )
                     FO(3) = CXZERO
                     FO(4) = DCMPLX( OMEGAP, RXZERO )
                  ELSE
                     RS2PA = RXONE / SQRT( RXTWO * PABS )
                     SPAZ = SQRT( PABS + P(3) )
                     FO(1) = OMEGAM * RS2PA / SPAZ
     &                     * DCMPLX( -P(1), -P(2) )
                     FO(2) = OMEGAM * RS2PA
     &                     * DCMPLX( SPAZ, RXZERO )
                     FO(3) = OMEGAP * RS2PA / SPAZ
     &                     * DCMPLX( -P(1), -P(2) )
                     FO(4) = OMEGAP * RS2PA
     &                     * DCMPLX( SPAZ, RXZERO )
                  END IF
               ELSE
                  IF ( PLAT .EQ. RXZERO ) THEN
                     FO(1) = DCMPLX( -OMEGAM, RXZERO )
                     FO(2) = CXZERO
                     FO(3) = DCMPLX( -OMEGAP, RXZERO )
                     FO(4) = CXZERO
                  ELSE
                     RS2PA = RXONE / SQRT( RXTWO * PABS )
                     SPAZ = SQRT( PABS - P(3) )
                     FO(1) = OMEGAM * RS2PA * SPAZ / PLAT
     &                     * DCMPLX( -P(1), -P(2) )
                     FO(2) = OMEGAM * RS2PA / SPAZ
     &                     * DCMPLX( PLAT, RXZERO )
                     FO(3) = OMEGAP * RS2PA * SPAZ / PLAT
     &                     * DCMPLX( -P(1), -P(2) )
                     FO(4) = OMEGAP * RS2PA / SPAZ
     &                     * DCMPLX( PLAT, RXZERO )
                  END IF
               END IF
            ELSE
               STOP 'OXXXXX:  FERMION HELICITY MUST BE +1,-1'
            END IF
         ELSE IF ( NSF .EQ. -1 ) THEN
            IF ( NHEL .EQ. 1 ) THEN
               IF ( P(3) .GE. RXZERO ) THEN
                  IF ( PLAT .EQ. RXZERO ) THEN
                     FO(1) = CXZERO
                     FO(2) = DCMPLX( OMEGAM, RXZERO )
                     FO(3) = CXZERO
                     FO(4) = DCMPLX( -OMEGAP, RXZERO )
                  ELSE
                     RS2PA = RXONE / SQRT( RXTWO * PABS )
                     SPAZ = SQRT( PABS + P(3) )
                     FO(1) = OMEGAM * RS2PA / SPAZ
     &                     * DCMPLX( -P(1), -P(2) )
                     FO(2) = OMEGAM * RS2PA
     &                     * DCMPLX( SPAZ, RXZERO )
                     FO(3) = -OMEGAP * RS2PA / SPAZ
     &                     * DCMPLX( -P(1), -P(2) )
                     FO(4) = -OMEGAP * RS2PA
     &                     * DCMPLX( SPAZ, RXZERO )
                  END IF
               ELSE
                  IF ( PLAT .EQ. RXZERO ) THEN
                     FO(1) = DCMPLX( -OMEGAM, RXZERO )
                     FO(2) = CXZERO
                     FO(3) = DCMPLX( OMEGAP, RXZERO )
                     FO(4) = CXZERO
                  ELSE
                     RS2PA = RXONE / SQRT( RXTWO * PABS )
                     SPAZ = SQRT( PABS - P(3) )
                     FO(1) = OMEGAM * RS2PA * SPAZ / PLAT
     &                     * DCMPLX( -P(1), -P(2) )
                     FO(2) = OMEGAM * RS2PA / SPAZ
     &                     * DCMPLX( PLAT, RXZERO )
                     FO(3) = -OMEGAP * RS2PA * SPAZ / PLAT
     &                     * DCMPLX( -P(1), -P(2) )
                     FO(4) = -OMEGAP * RS2PA / SPAZ
     &                     * DCMPLX( PLAT, RXZERO )
                  END IF
               END IF
            ELSE IF ( NHEL .EQ. -1 ) THEN
               IF ( P(3) .GE. RXZERO ) THEN
                  IF ( PLAT .EQ. RXZERO ) THEN
                     FO(1) = DCMPLX( -OMEGAP, RXZERO )
                     FO(2) = CXZERO
                     FO(3) = DCMPLX( OMEGAM, RXZERO )
                     FO(4) = CXZERO
                  ELSE
                     RS2PA = RXONE / SQRT( RXTWO * PABS )
                     SPAZ = SQRT( PABS + P(3) )
                     FO(1) = -OMEGAP * RS2PA
     &                     * DCMPLX( SPAZ, RXZERO )
                     FO(2) = -OMEGAP * RS2PA / SPAZ
     &                     * DCMPLX( P(1), -P(2) )
                     FO(3) = OMEGAM * RS2PA
     &                     * DCMPLX( SPAZ, RXZERO )
                     FO(4) = OMEGAM * RS2PA / SPAZ
     &                     * DCMPLX( P(1), -P(2) )
                  END IF
               ELSE
                  IF ( PLAT .EQ. RXZERO ) THEN
                     FO(1) = CXZERO
                     FO(2) = DCMPLX( -OMEGAP, RXZERO )
                     FO(3) = CXZERO
                     FO(4) = DCMPLX( OMEGAM, RXZERO )
                  ELSE
                     RS2PA = RXONE / SQRT( RXTWO * PABS )
                     SPAZ = SQRT( PABS - P(3) )
                     FO(1) = -OMEGAP * RS2PA / SPAZ
     &                     * DCMPLX( PLAT, RXZERO )
                     FO(2) = -OMEGAP * RS2PA * SPAZ / PLAT
     &                     * DCMPLX( P(1), -P(2) )
                     FO(3) = OMEGAM * RS2PA / SPAZ
     &                     * DCMPLX( PLAT, RXZERO )
                     FO(4) = OMEGAM * RS2PA * SPAZ / PLAT
     &                     * DCMPLX( P(1), -P(2) )
                  END IF
               END IF
            ELSE
               STOP 'OXXXXX:  FERMION HELICITY MUST BE +1,-1'
            END IF
         ELSE
            STOP 'OXXXXX:  FERMION TYPE MUST BE +1,-1'
         END IF
C
C do massless case
C
      ELSE
         IF ( NSF .EQ. 1 ) THEN
            IF ( NHEL .EQ. 1 ) THEN
               IF ( P(3) .GE. RXZERO ) THEN
                  IF ( PLAT .EQ. RXZERO ) THEN
                     FO(1) = DCMPLX( OMEGAP, RXZERO )
                     FO(2) = CXZERO
                     FO(3) = CXZERO
                     FO(4) = CXZERO
                  ELSE
                     SPAZ = SQRT( PABS + P(3) )
                     FO(1) = DCMPLX( SPAZ, RXZERO )
                     FO(2) = RXONE / SPAZ
     &                     * DCMPLX( P(1), -P(2) )
                     FO(3) = CXZERO
                     FO(4) = CXZERO
                  END IF
               ELSE
                  IF ( PLAT .EQ. RXZERO ) THEN
                     FO(1) = CXZERO
                     FO(2) = DCMPLX( OMEGAP, RXZERO )
                     FO(3) = CXZERO
                     FO(4) = CXZERO
                  ELSE
                     SPAZ = SQRT( PABS - P(3) )
                     FO(1) = RXONE / SPAZ
     &                     * DCMPLX( PLAT, RXZERO )
                     FO(2) = SPAZ / PLAT
     &                     * DCMPLX( P(1), -P(2) )
                     FO(3) = CXZERO
                     FO(4) = CXZERO
                  END IF
               END IF
            ELSE IF ( NHEL .EQ. -1 ) THEN
               IF ( P(3) .GE. RXZERO ) THEN
                  IF ( PLAT .EQ. RXZERO ) THEN
                     FO(1) = CXZERO
                     FO(2) = CXZERO
                     FO(3) = CXZERO
                     FO(4) = DCMPLX( OMEGAP, RXZERO )
                  ELSE
                     SPAZ = SQRT( PABS + P(3) )
                     FO(1) = CXZERO
                     FO(2) = CXZERO
                     FO(3) = RXONE / SPAZ
     &                     * DCMPLX( -P(1), -P(2) )
                     FO(4) = DCMPLX( SPAZ, RXZERO )
                  END IF
               ELSE
                  IF ( PLAT .EQ. RXZERO ) THEN
                     FO(1) = CXZERO
                     FO(2) = CXZERO
                     FO(3) = DCMPLX( -OMEGAP, RXZERO )
                     FO(4) = CXZERO
                  ELSE
                     SPAZ = SQRT( PABS - P(3) )
                     FO(1) = CXZERO
                     FO(2) = CXZERO
                     FO(3) = SPAZ / PLAT
     &                     * DCMPLX( -P(1), -P(2) )
                     FO(4) = RXONE / SPAZ
     &                     * DCMPLX( PLAT, RXZERO )
                  END IF
               END IF
            ELSE
               STOP 'OXXXXX:  FERMION HELICITY MUST BE +1,-1'
            END IF
         ELSE IF ( NSF .EQ. -1 ) THEN
            IF ( NHEL .EQ. 1 ) THEN
               IF ( P(3) .GE. RXZERO ) THEN
                  IF ( PLAT .EQ. RXZERO ) THEN
                     FO(1) = CXZERO
                     FO(2) = CXZERO
                     FO(3) = CXZERO
                     FO(4) = DCMPLX( -OMEGAP, RXZERO )
                  ELSE
                     SPAZ = SQRT( PABS + P(3) )
                     FO(1) = CXZERO
                     FO(2) = CXZERO
                     FO(3) = -RXONE / SPAZ
     &                     * DCMPLX( -P(1), -P(2) )
                     FO(4) = DCMPLX( -SPAZ, RXZERO )
                  END IF
               ELSE
                  IF ( PLAT .EQ. RXZERO ) THEN
                     FO(1) = CXZERO
                     FO(2) = CXZERO
                     FO(3) = DCMPLX( OMEGAP, RXZERO )
                     FO(4) = CXZERO
                  ELSE
                     SPAZ = SQRT( PABS - P(3) )
                     FO(1) = CXZERO
                     FO(2) = CXZERO
                     FO(3) = -SPAZ / PLAT
     &                     * DCMPLX( -P(1), -P(2) )
                     FO(4) = -RXONE / SPAZ
     &                     * DCMPLX( PLAT, RXZERO )
                  END IF
               END IF
            ELSE IF ( NHEL .EQ. -1 ) THEN
               IF ( P(3) .GE. RXZERO ) THEN
                  IF ( PLAT .EQ. RXZERO ) THEN
                     FO(1) = DCMPLX( -OMEGAP, RXZERO )
                     FO(2) = CXZERO
                     FO(3) = CXZERO
                     FO(4) = CXZERO
                  ELSE
                     SPAZ = SQRT( PABS + P(3) )
                     FO(1) = DCMPLX( -SPAZ, RXZERO )
                     FO(2) = -RXONE / SPAZ
     &                     * DCMPLX( P(1), -P(2) )
                     FO(3) = CXZERO
                     FO(4) = CXZERO
                  END IF
               ELSE
                  IF ( PLAT .EQ. RXZERO ) THEN
                     FO(1) = CXZERO
                     FO(2) = DCMPLX( -OMEGAP, RXZERO )
                     FO(3) = CXZERO
                     FO(4) = CXZERO
                  ELSE
                     SPAZ = SQRT( PABS - P(3) )
                     FO(1) = -RXONE / SPAZ
     &                     * DCMPLX( PLAT, RXZERO )
                     FO(2) = -SPAZ / PLAT
     &                     * DCMPLX( P(1), -P(2) )
                     FO(3) = CXZERO
                     FO(4) = CXZERO
                  END IF
               END IF
            ELSE
               STOP 'OXXXXX:  FERMION HELICITY MUST BE +1,-1'
            END IF
         ELSE
            STOP 'OXXXXX:  FERMION TYPE MUST BE +1,-1'
         END IF
      END IF
C
C done
C
      RETURN
      END
C
C ----------------------------------------------------------------------
C
      SUBROUTINE ROTXXX(P,Q , PROT)
C
C this subroutine performs the spacial rotation of a four-momentum.     
C the momentum p is assumed to be given in the frame where the spacial  
C component of q points the positive z-axis.  prot is the momentum p    
C rotated to the frame where q is given.                                
C                                                                       
C input:                                                                
C       real    p(0:3)         : four-momentum p in q(1)=q(2)=0 frame   
C       real    q(0:3)         : four-momentum q in the rotated frame   
C                                                                       
C output:                                                               
C       real    prot(0:3)      : four-momentum p in the rotated frame   
C
#if defined(CERNLIB_IMPNONE)
      IMPLICIT NONE
#endif
      REAL*8    P(0:3),Q(0:3),PROT(0:3),QT2,QT,PSGN,QQ,P1
C
      REAL*8 RXZERO, RXONE
      PARAMETER( RXZERO=0.0D0, RXONE=1.0D0 )
C
      PROT(0) = P(0)
C
      QT2=Q(1)**2+Q(2)**2
C
      IF ( QT2 .EQ. RXZERO ) THEN
          IF ( Q(3) .EQ. RXZERO ) THEN
             PROT(1) = P(1)
             PROT(2) = P(2)
             PROT(3) = P(3)
          ELSE
             PSGN=DSIGN(RXONE,Q(3))
             PROT(1) = P(1)*PSGN
             PROT(2) = P(2)*PSGN
             PROT(3) = P(3)*PSGN
          ENDIF
      ELSE
          QQ=SQRT(QT2+Q(3)**2)
          QT=SQRT(QT2)
          P1=P(1)
          PROT(1) = Q(1)*Q(3)/QQ/QT*P1 -Q(2)/QT*P(2) +Q(1)/QQ*P(3)
          PROT(2) = Q(2)*Q(3)/QQ/QT*P1 +Q(1)/QT*P(2) +Q(2)/QQ*P(3)
          PROT(3) =          -QT/QQ*P1               +Q(3)/QQ*P(3)
      ENDIF
C
      RETURN
      END
C ======================================================================
C
      SUBROUTINE SSSSXX(S1,S2,S3,S4,G , VERTEX)
C
C This subroutine computes an amplitude of the four-scalar coupling.    
C                                                                       
C INPUT:                                                                
C       complex S1(3)          : first  scalar                        S1
C       complex S2(3)          : second scalar                        S2
C       complex S3(3)          : third  scalar                        S3
C       complex S4(3)          : fourth scalar                        S4
C       real    G              : coupling constant                 GHHHH
C                                                                       
C OUTPUT:                                                               
C       complex VERTEX         : amplitude            Gamma(S1,S2,S3,S4)
C
#if defined(CERNLIB_IMPNONE)
      IMPLICIT NONE
#endif
      COMPLEX*16 S1(3),S2(3),S3(3),S4(3),VERTEX
      REAL*8     G
C
      VERTEX = G*S1(1)*S2(1)*S3(1)*S4(1)
C
      RETURN
      END
C
C ======================================================================
C
      SUBROUTINE SSSXXX(S1,S2,S3,G , VERTEX)
C
C This subroutine computes an amplitude of the three-scalar coupling.   
C                                                                       
C INPUT:                                                                
C       complex S1(3)          : first  scalar                        S1
C       complex S2(3)          : second scalar                        S2
C       complex S3(3)          : third  scalar                        S3
C       real    G              : coupling constant                  GHHH
C                                                                       
C OUTPUT:                                                               
C       complex VERTEX         : amplitude               Gamma(S1,S2,S3)
C
#if defined(CERNLIB_IMPNONE)
      IMPLICIT NONE
#endif
      COMPLEX*16 S1(3),S2(3),S3(3),VERTEX
      REAL*8    G
C
      VERTEX = G*S1(1)*S2(1)*S3(1)
C
      RETURN
      END
C
C
C ----------------------------------------------------------------------
C
      SUBROUTINE SXXXXX(P,NSS , SC)
C
C This subroutine computes a complex SCALAR wavefunction.               
C                                                                       
C INPUT:                                                                
C       real    P(0:3)         : four-momentum of scalar boson          
C       integer NSS  = -1 or 1 : +1 for final, -1 for initial           
C                                                                       
C OUTPUT:                                                               
C       complex SC(3)          : scalar wavefunction                   S
C
#if defined(CERNLIB_IMPNONE)
      IMPLICIT NONE
#endif
      COMPLEX*16 SC(3)
      REAL*8    P(0:3)
      INTEGER NSS
C
      SC(1) = DCMPLX( 1.0 )
      SC(2) = DCMPLX(P(0),P(3))*NSS
      SC(3) = DCMPLX(P(1),P(2))*NSS
C
      RETURN
      END
C
C ======================================================================
C
      SUBROUTINE VSSXXX(VC,S1,S2,G , VERTEX)
C
C this subroutine computes an amplitude from the vector-scalar-scalar   
C coupling.  the coupling is absent in the minimal sm in unitary gauge. 
C                                                                       
C       complex vc(6)          : input  vector                        v 
C       complex s1(3)          : first  scalar                        s1
C       complex s2(3)          : second scalar                        s2
C       complex g              : coupling constant (s1 charge)          
C                                                                       
C examples of the coupling constant g for susy particles are as follows:
C   -----------------------------------------------------------         
C   |    s1    | (q,i3) of s1  ||   v=a   |   v=z   |   v=w   |         
C   -----------------------------------------------------------         
C   | nu~_l    | (  0  , +1/2) ||   ---   |  gzn(1) |  gwf(1) |         
C   | e~_l     | ( -1  , -1/2) ||  gal(1) |  gzl(1) |  gwf(1) |         
C   | u~_l     | (+2/3 , +1/2) ||  gau(1) |  gzu(1) |  gwf(1) |         
C   | d~_l     | (-1/3 , -1/2) ||  gad(1) |  gzd(1) |  gwf(1) |         
C   -----------------------------------------------------------         
C   | e~_r-bar | ( +1  ,  0  ) || -gal(2) | -gzl(2) | -gwf(2) |         
C   | u~_r-bar | (-2/3 ,  0  ) || -gau(2) | -gzu(2) | -gwf(2) |         
C   | d~_r-bar | (+1/3 ,  0  ) || -gad(2) | -gzd(2) | -gwf(2) |         
C   -----------------------------------------------------------         
C where the s1 charge is defined by the flowing-out quantum number.     
C                                                                       
C output:                                                               
C       complex vertex         : amplitude                gamma(v,s1,s2)
C
#if defined(CERNLIB_IMPNONE)
      IMPLICIT NONE
#endif
      COMPLEX*16 VC(6),S1(3),S2(3),VERTEX,G
      REAL*8    P(0:3)
C
      P(0)=DBLE( S1(2)-S2(2))
      P(1)=DBLE( S1(3)-S2(3))
      P(2)=DIMAG(S1(3)-S2(3))
      P(3)=DIMAG(S1(2)-S2(2))
C
      VERTEX = G*S1(1)*S2(1)
     &        *(VC(1)*P(0)-VC(2)*P(1)-VC(3)*P(2)-VC(4)*P(3))
C
      RETURN
      END
C
      SUBROUTINE VVSSXX(V1,V2,S1,S2,G , VERTEX)
C
C This subroutine computes an amplitude of the vector-vector-scalar-    
C scalar coupling.                                                      
C                                                                       
C INPUT:                                                                
C       complex V1(6)          : first  vector                        V1
C       complex V2(6)          : second vector                        V2
C       complex S1(3)          : first  scalar                        S1
C       complex S2(3)          : second scalar                        S2
C       real    G              : coupling constant                 GVVHH
C                                                                       
C OUTPUT:                                                               
C       complex VERTEX         : amplitude            Gamma(V1,V2,S1,S2)
C
#if defined(CERNLIB_IMPNONE)
      IMPLICIT NONE
#endif
      COMPLEX*16 V1(6),V2(6),S1(3),S2(3),VERTEX
      REAL*8    G
C
      VERTEX = G*S1(1)*S2(1)
     &        *(V1(1)*V2(1)-V1(2)*V2(2)-V1(3)*V2(3)-V1(4)*V2(4))
C
      RETURN
      END
C
C
C ======================================================================
C
      SUBROUTINE VVSXXX(V1,V2,SC,G , VERTEX)
C
C this subroutine computes an amplitude of the vector-vector-scalar     
C coupling.                                                             
C                                                                       
C input:                                                                
C       complex v1(6)          : first  vector                        v1
C       complex v2(6)          : second vector                        v2
C       complex sc(3)          : input  scalar                        s 
C       real    g              : coupling constant                  gvvh
C                                                                       
C output:                                                               
C       complex vertex         : amplitude                gamma(v1,v2,s)
C
#if defined(CERNLIB_IMPNONE)
      IMPLICIT NONE
#endif
      COMPLEX*16 V1(6),V2(6),SC(3),VERTEX
      REAL*8    G
C
      VERTEX = G*SC(1)*(V1(1)*V2(1)-V1(2)*V2(2)-V1(3)*V2(3)-V1(4)*V2(4))
C
      RETURN
      END
C
C ======================================================================
C
      SUBROUTINE VVVXXX(WM,WP,W3,G , VERTEX)
C
C this subroutine computes an amplitude of the three-point coupling of  
C the gauge bosons.                                                     
C                                                                       
C input:                                                                
C       complex wm(6)          : vector               flow-out w-       
C       complex wp(6)          : vector               flow-out w+       
C       complex w3(6)          : vector               j3 or a    or z   
C       real    g              : coupling constant    gw or gwwa or gwwz
C                                                                       
C output:                                                               
C       complex vertex         : amplitude               gamma(wm,wp,w3)
C
#if defined(CERNLIB_IMPNONE)
      IMPLICIT NONE
#endif
      COMPLEX*16 WM(6),WP(6),W3(6),VERTEX, 
     &        XV1,XV2,XV3,V12,V23,V31,P12,P13,P21,P23,P31,P32
      REAL*8    PWM(0:3),PWP(0:3),PW3(0:3),G
C
      REAL*8 RXZERO, RTENTH
      PARAMETER( RXZERO=0.0D0, RTENTH=0.1D0 )
C
      PWM(0)=DBLE( WM(5))
      PWM(1)=DBLE( WM(6))
      PWM(2)=DIMAG(WM(6))
      PWM(3)=DIMAG(WM(5))
      PWP(0)=DBLE( WP(5))
      PWP(1)=DBLE( WP(6))
      PWP(2)=DIMAG(WP(6))
      PWP(3)=DIMAG(WP(5))
      PW3(0)=DBLE( W3(5))
      PW3(1)=DBLE( W3(6))
      PW3(2)=DIMAG(W3(6))
      PW3(3)=DIMAG(W3(5))
C
      V12=WM(1)*WP(1)-WM(2)*WP(2)-WM(3)*WP(3)-WM(4)*WP(4)
      V23=WP(1)*W3(1)-WP(2)*W3(2)-WP(3)*W3(3)-WP(4)*W3(4)
      V31=W3(1)*WM(1)-W3(2)*WM(2)-W3(3)*WM(3)-W3(4)*WM(4)
      XV1=RXZERO
      XV2=RXZERO
      XV3=RXZERO
      IF ( ABS(WM(1)) .NE. RXZERO ) THEN
         IF (ABS(WM(1)).GE.MAX(ABS(WM(2)),ABS(WM(3)),ABS(WM(4)))
     $        *RTENTH)
     &      XV1=PWM(0)/WM(1)
      ENDIF
      IF ( ABS(WP(1)) .NE. RXZERO) THEN
         IF (ABS(WP(1)).GE.MAX(ABS(WP(2)),ABS(WP(3)),ABS(WP(4)))
     $        *RTENTH)
     &      XV2=PWP(0)/WP(1)
      ENDIF
      IF ( ABS(W3(1)) .NE. RXZERO) THEN
         IF ( ABS(W3(1)).GE.MAX(ABS(W3(2)),ABS(W3(3)),ABS(W3(4)))
     $        *RTENTH)
     &      XV3=PW3(0)/W3(1)
      ENDIF
      P12= (PWM(0)-XV1*WM(1))*WP(1)-(PWM(1)-XV1*WM(2))*WP(2)
     &    -(PWM(2)-XV1*WM(3))*WP(3)-(PWM(3)-XV1*WM(4))*WP(4)
      P13= (PWM(0)-XV1*WM(1))*W3(1)-(PWM(1)-XV1*WM(2))*W3(2)
     &    -(PWM(2)-XV1*WM(3))*W3(3)-(PWM(3)-XV1*WM(4))*W3(4)
      P21= (PWP(0)-XV2*WP(1))*WM(1)-(PWP(1)-XV2*WP(2))*WM(2)
     &    -(PWP(2)-XV2*WP(3))*WM(3)-(PWP(3)-XV2*WP(4))*WM(4)
      P23= (PWP(0)-XV2*WP(1))*W3(1)-(PWP(1)-XV2*WP(2))*W3(2)
     &    -(PWP(2)-XV2*WP(3))*W3(3)-(PWP(3)-XV2*WP(4))*W3(4)
      P31= (PW3(0)-XV3*W3(1))*WM(1)-(PW3(1)-XV3*W3(2))*WM(2)
     &    -(PW3(2)-XV3*W3(3))*WM(3)-(PW3(3)-XV3*W3(4))*WM(4)
      P32= (PW3(0)-XV3*W3(1))*WP(1)-(PW3(1)-XV3*W3(2))*WP(2)
     &    -(PW3(2)-XV3*W3(3))*WP(3)-(PW3(3)-XV3*W3(4))*WP(4)
C
      VERTEX = -(V12*(P13-P23)+V23*(P21-P31)+V31*(P32-P12))*G
C
      RETURN
      END
C
C
C	Subroutine returns the value of evaluated
C	helicity basis boson polarisation wavefunction.
C	Replaces the HELAS routine VXXXXX
C
C	Adam Duff,  1992 September 3
C	<duff@phenom.physics.wisc.edu>
C
      SUBROUTINE VXXXXX(P,VMASS,NHEL,NSV,VC)
C
C          P            IN: BOSON FOUR MOMENTUM
C          VMASS        IN: BOSON MASS
C          NHEL         IN: BOSON HELICITY
C          NSV          IN: INCOMING (-1) OR OUTGOING (+1)
C          VC           OUT: BOSON WAVEFUNCTION
C
C declare input/output variables
C
#if defined(CERNLIB_IMPNONE)
      IMPLICIT NONE
#endif
      COMPLEX*16 VC(6)
      INTEGER*4 NHEL, NSV
      REAL*8 P(0:3), VMASS
C
C declare local variables
C
      REAL*8 RXZERO, RXONE, RXTWO
      PARAMETER( RXZERO=0.0D0, RXONE=1.0D0, RXTWO=2.0D0 )
      REAL*8 PLAT, PABS, RS2, RPLAT, RPABS, RDEN
      COMPLEX*16 CXZERO
      LOGICAL FIRST
      SAVE CXZERO,FIRST
      DATA FIRST/.TRUE./
C
C          Fix compilation with g77
      IF(FIRST) THEN
        FIRST=.FALSE.
        CXZERO=DCMPLX( RXZERO, RXZERO )
      ENDIF
C
C define internal/external momenta
C
      IF ( NSV**2 .NE. 1 ) THEN
         STOP 'VXXXXX:  NSV IS NOT ONE OF -1, +1'
      END IF
C
      RS2 = SQRT( RXONE / RXTWO )
      VC(5) = DCMPLX( P(0), P(3) ) * NSV
      VC(6) = DCMPLX( P(1), P(2) ) * NSV
      PLAT = SQRT( P(1)**2 + P(2)**2 )
      PABS = SQRT( P(1)**2 + P(2)**2 + P(3)**2 )
C
C calculate polarisation four vectors
C
      IF ( NHEL**2 .EQ. 1 ) THEN
         IF ( (PABS .EQ. RXZERO) .OR. (PLAT .EQ. RXZERO) ) THEN
            VC(1) = CXZERO
            VC(2) = DCMPLX( -NHEL * RS2 * DSIGN( RXONE, P(3) ), RXZERO )
            VC(3) = DCMPLX( RXZERO, NSV * RS2 )
            VC(4) = CXZERO
         ELSE
            RPLAT = RXONE / PLAT
            RPABS = RXONE / PABS
            VC(1) = CXZERO
            VC(2) = DCMPLX( -NHEL * RS2 * RPABS * RPLAT * P(1) * P(3),
     &                     -NSV * RS2 * RPLAT * P(2) )
            VC(3) = DCMPLX( -NHEL * RS2 * RPABS * RPLAT * P(2) * P(3),
     &                     NSV * RS2 * RPLAT * P(1) )
            VC(4) = DCMPLX( NHEL * RS2 * RPABS * PLAT,
     &                     RXZERO )
         END IF
      ELSE IF ( NHEL .EQ. 0 ) THEN
         IF ( VMASS .GT. RXZERO ) THEN
            IF ( PABS .EQ. RXZERO ) THEN
               VC(1) = CXZERO
               VC(2) = CXZERO
               VC(3) = CXZERO
               VC(4) = DCMPLX( RXONE, RXZERO )
            ELSE
               RDEN = P(0) / ( VMASS * PABS )
               VC(1) = DCMPLX( PABS / VMASS, RXZERO )
               VC(2) = DCMPLX( RDEN * P(1), RXZERO )
               VC(3) = DCMPLX( RDEN * P(2), RXZERO )
               VC(4) = DCMPLX( RDEN * P(3), RXZERO )
            END IF
         ELSE
            STOP  'VXXXXX: NHEL = 0 IS ONLY FOR MASSIVE BOSONS'
         END IF
      ELSE IF ( NHEL .EQ. 4 ) THEN
         IF ( VMASS .GT. RXZERO ) THEN
            RDEN = RXONE / VMASS
            VC(1) = DCMPLX( RDEN * P(0), RXZERO )
            VC(2) = DCMPLX( RDEN * P(1), RXZERO )
            VC(3) = DCMPLX( RDEN * P(2), RXZERO )
            VC(4) = DCMPLX( RDEN * P(3), RXZERO )
         ELSEIF (VMASS .EQ. RXZERO) THEN
            RDEN = RXONE / P(0)
            VC(1) = DCMPLX( RDEN * P(0), RXZERO )
            VC(2) = DCMPLX( RDEN * P(1), RXZERO )
            VC(3) = DCMPLX( RDEN * P(2), RXZERO )
            VC(4) = DCMPLX( RDEN * P(3), RXZERO )
         ELSE
            STOP 'VXXXXX: NHEL = 4 IS ONLY FOR M>=0'
         END IF
      ELSE
         STOP 'VXXXXX:  NHEL IS NOT ONE OF -1, 0, 1 OR 4'
      END IF
C
C done
C
      RETURN
      END
C
C ----------------------------------------------------------------------
C
      SUBROUTINE W3W3XX(WM,W31,WP,W32,G31,G32,WMASS,WWIDTH , VERTEX)
C
C this subroutine computes an amplitude of the four-point coupling of   
C the w-, w+ and two w3/z/a.  the amplitude includes the contributions  
C of w exchange diagrams.  the internal w propagator is given in unitary
C gauge.  if one sets wmass=0.0, then the gggg vertex is given (see sect
C 2.9.1 of the manual).
C                                                                       
C input:                                                                
C       complex wm(0:3)        : flow-out w-                         wm 
C       complex w31(0:3)       : first    w3/z/a                     w31
C       complex wp(0:3)        : flow-out w+                         wp 
C       complex w32(0:3)       : second   w3/z/a                     w32
C       real    g31            : coupling of w31 with w-/w+             
C       real    g32            : coupling of w32 with w-/w+             
C                                                  (see the table below)
C       real    wmass          : mass  of w                             
C       real    wwidth         : width of w                             
C                                                                       
C the possible sets of the inputs are as follows:                       
C   -------------------------------------------                         
C   |  wm  |  w31 |  wp  |  w32 |  g31 |  g32 |                         
C   -------------------------------------------                         
C   |  w-  |  w3  |  w+  |  w3  |  gw  |  gw  |                         
C   |  w-  |  w3  |  w+  |  z   |  gw  | gwwz |                         
C   |  w-  |  w3  |  w+  |  a   |  gw  | gwwa |                         
C   |  w-  |  z   |  w+  |  z   | gwwz | gwwz |                         
C   |  w-  |  z   |  w+  |  a   | gwwz | gwwa |                         
C   |  w-  |  a   |  w+  |  a   | gwwa | gwwa |                         
C   -------------------------------------------                         
C where all the bosons are defined by the flowing-out quantum number.   
C                                                                       
C output:                                                               
C       complex vertex         : amplitude          gamma(wm,w31,wp,w32)
C
#if defined(CERNLIB_IMPNONE)
      IMPLICIT NONE
#endif
      COMPLEX*16    WM(6),W31(6),WP(6),W32(6),VERTEX
      COMPLEX*16 DV1(0:3),DV2(0:3),DV3(0:3),DV4(0:3),DVERTX,
     &           V12,V13,V14,V23,V24,V34
      REAL*8     G31,G32,WMASS,WWIDTH
C
      REAL*8 RXZERO, RXONE
      PARAMETER( RXZERO=0.0D0, RXONE=1.0D0 )

      DV1(0)=DCMPLX(WM(1))
      DV1(1)=DCMPLX(WM(2))
      DV1(2)=DCMPLX(WM(3))
      DV1(3)=DCMPLX(WM(4))
      DV2(0)=DCMPLX(W31(1))
      DV2(1)=DCMPLX(W31(2))
      DV2(2)=DCMPLX(W31(3))
      DV2(3)=DCMPLX(W31(4))
      DV3(0)=DCMPLX(WP(1))
      DV3(1)=DCMPLX(WP(2))
      DV3(2)=DCMPLX(WP(3))
      DV3(3)=DCMPLX(WP(4))
      DV4(0)=DCMPLX(W32(1))
      DV4(1)=DCMPLX(W32(2))
      DV4(2)=DCMPLX(W32(3))
      DV4(3)=DCMPLX(W32(4))
C
      IF ( DBLE(WMASS) .NE. RXZERO ) THEN
C         dm2inv = r_one / dmw2
C
         V12= DV1(0)*DV2(0)-DV1(1)*DV2(1)-DV1(2)*DV2(2)-DV1(3)*DV2(3)
         V13= DV1(0)*DV3(0)-DV1(1)*DV3(1)-DV1(2)*DV3(2)-DV1(3)*DV3(3)
         V14= DV1(0)*DV4(0)-DV1(1)*DV4(1)-DV1(2)*DV4(2)-DV1(3)*DV4(3)
         V23= DV2(0)*DV3(0)-DV2(1)*DV3(1)-DV2(2)*DV3(2)-DV2(3)*DV3(3)
         V24= DV2(0)*DV4(0)-DV2(1)*DV4(1)-DV2(2)*DV4(2)-DV2(3)*DV4(3)
         V34= DV3(0)*DV4(0)-DV3(1)*DV4(1)-DV3(2)*DV4(2)-DV3(3)*DV4(3)
C
         DVERTX = V12*V34 +V14*V23 -2.D0*V13*V24
C
         VERTEX = DCMPLX( DVERTX ) * (G31*G32)
C
      ELSE
         V12= DV1(0)*DV2(0)-DV1(1)*DV2(1)-DV1(2)*DV2(2)-DV1(3)*DV2(3)
         V13= DV1(0)*DV3(0)-DV1(1)*DV3(1)-DV1(2)*DV3(2)-DV1(3)*DV3(3)
         V14= DV1(0)*DV4(0)-DV1(1)*DV4(1)-DV1(2)*DV4(2)-DV1(3)*DV4(3)
         V23= DV2(0)*DV3(0)-DV2(1)*DV3(1)-DV2(2)*DV3(2)-DV2(3)*DV3(3)
         V24= DV2(0)*DV4(0)-DV2(1)*DV4(1)-DV2(2)*DV4(2)-DV2(3)*DV4(3)
         V34= DV3(0)*DV4(0)-DV3(1)*DV4(1)-DV3(2)*DV4(2)-DV3(3)*DV4(3)
C

         DVERTX = V14*V23 -V13*V24
C
         VERTEX = DCMPLX( DVERTX ) * (G31*G32)
      END IF
C
      RETURN
      END
C
C ======================================================================
C
      SUBROUTINE WWWWXX(WM1,WP1,WM2,WP2,GWWA,GWWZ,ZMASS,ZWIDTH , VERTEX)
C
C this subroutine computes an amplitude of the four-point w-/w+         
C coupling, including the contributions of photon and z exchanges.  the 
C photon propagator is given in feynman gauge and the z propagator is   
C given in unitary gauge.                                               
C                                                                       
C input:                                                                
C       complex wm1(0:3)       : first  flow-out w-                  wm1
C       complex wp1(0:3)       : first  flow-out w+                  wp1
C       complex wm2(0:3)       : second flow-out w-                  wm2
C       complex wp2(0:3)       : second flow-out w+                  wp2
C       real    gwwa           : coupling constant of w and a       gwwa
C       real    gwwz           : coupling constant of w and z       gwwz
C       real    zmass          : mass  of z                             
C       real    zwidth         : width of z                             
C                                                                       
C output:                                                               
C       complex vertex         : amplitude        gamma(wm1,wp1,wm2,wp2)
C
#if defined(CERNLIB_IMPNONE)
      IMPLICIT NONE
#endif
      COMPLEX*16    WM1(6),WP1(6),WM2(6),WP2(6),VERTEX
      COMPLEX*16 DV1(0:3),DV2(0:3),DV3(0:3),DV4(0:3),
     &           J12(0:3),J34(0:3),J14(0:3),J32(0:3),DVERTX,
     &           SV1,SV2,SV3,SV4,TV1,TV2,TV3,TV4,DZS,DZT,
     &           V12,V13,V14,V23,V24,V34,JS12,JS34,JS14,JS32,JS,JT
      REAL*8       PWM1(0:3),PWP1(0:3),PWM2(0:3),PWP2(0:3),
     &           GWWA,GWWZ,ZMASS,ZWIDTH
      REAL*8     Q(0:3),K(0:3),DP1(0:3),DP2(0:3),DP3(0:3),DP4(0:3),
     &           DGWWA2,DGWWZ2,DGW2,DMZ,DWIDTH,S,T,DAS,DAT
C
      REAL*8 RXZERO, RXONE, RXTWO
      PARAMETER( RXZERO=0.0D0, RXONE=1.0D0, RXTWO=2.0D0 )
C
      PWM1(0)=DBLE( WM1(5))
      PWM1(1)=DBLE( WM1(6))
      PWM1(2)=DIMAG(WM1(6))
      PWM1(3)=DIMAG(WM1(5))
      PWP1(0)=DBLE( WP1(5))
      PWP1(1)=DBLE( WP1(6))
      PWP1(2)=DIMAG(WP1(6))
      PWP1(3)=DIMAG(WP1(5))
      PWM2(0)=DBLE( WM2(5))
      PWM2(1)=DBLE( WM2(6))
      PWM2(2)=DIMAG(WM2(6))
      PWM2(3)=DIMAG(WM2(5))
      PWP2(0)=DBLE( WP2(5))
      PWP2(1)=DBLE( WP2(6))
      PWP2(2)=DIMAG(WP2(6))
      PWP2(3)=DIMAG(WP2(5))
C
      DV1(0)=DCMPLX(WM1(1))
      DV1(1)=DCMPLX(WM1(2))
      DV1(2)=DCMPLX(WM1(3))
      DV1(3)=DCMPLX(WM1(4))
      DP1(0)=DBLE(PWM1(0))
      DP1(1)=DBLE(PWM1(1))
      DP1(2)=DBLE(PWM1(2))
      DP1(3)=DBLE(PWM1(3))
      DV2(0)=DCMPLX(WP1(1))
      DV2(1)=DCMPLX(WP1(2))
      DV2(2)=DCMPLX(WP1(3))
      DV2(3)=DCMPLX(WP1(4))
      DP2(0)=DBLE(PWP1(0))
      DP2(1)=DBLE(PWP1(1))
      DP2(2)=DBLE(PWP1(2))
      DP2(3)=DBLE(PWP1(3))
      DV3(0)=DCMPLX(WM2(1))
      DV3(1)=DCMPLX(WM2(2))
      DV3(2)=DCMPLX(WM2(3))
      DV3(3)=DCMPLX(WM2(4))
      DP3(0)=DBLE(PWM2(0))
      DP3(1)=DBLE(PWM2(1))
      DP3(2)=DBLE(PWM2(2))
      DP3(3)=DBLE(PWM2(3))
      DV4(0)=DCMPLX(WP2(1))
      DV4(1)=DCMPLX(WP2(2))
      DV4(2)=DCMPLX(WP2(3))
      DV4(3)=DCMPLX(WP2(4))
      DP4(0)=DBLE(PWP2(0))
      DP4(1)=DBLE(PWP2(1))
      DP4(2)=DBLE(PWP2(2))
      DP4(3)=DBLE(PWP2(3))
      DGWWA2=DBLE(GWWA)**2
      DGWWZ2=DBLE(GWWZ)**2
      DGW2  =DGWWA2+DGWWZ2
      DMZ   =DBLE(ZMASS)
      DWIDTH=DBLE(ZWIDTH)
C
      V12= DV1(0)*DV2(0)-DV1(1)*DV2(1)-DV1(2)*DV2(2)-DV1(3)*DV2(3)
      V13= DV1(0)*DV3(0)-DV1(1)*DV3(1)-DV1(2)*DV3(2)-DV1(3)*DV3(3)
      V14= DV1(0)*DV4(0)-DV1(1)*DV4(1)-DV1(2)*DV4(2)-DV1(3)*DV4(3)
      V23= DV2(0)*DV3(0)-DV2(1)*DV3(1)-DV2(2)*DV3(2)-DV2(3)*DV3(3)
      V24= DV2(0)*DV4(0)-DV2(1)*DV4(1)-DV2(2)*DV4(2)-DV2(3)*DV4(3)
      V34= DV3(0)*DV4(0)-DV3(1)*DV4(1)-DV3(2)*DV4(2)-DV3(3)*DV4(3)
C
      Q(0)=DP1(0)+DP2(0)
      Q(1)=DP1(1)+DP2(1)
      Q(2)=DP1(2)+DP2(2)
      Q(3)=DP1(3)+DP2(3)
      K(0)=DP1(0)+DP4(0)
      K(1)=DP1(1)+DP4(1)
      K(2)=DP1(2)+DP4(2)
      K(3)=DP1(3)+DP4(3)
C
      S=Q(0)**2-Q(1)**2-Q(2)**2-Q(3)**2
      T=K(0)**2-K(1)**2-K(2)**2-K(3)**2
C
      DAS=-RXONE/S
      DAT=-RXONE/T
      DZS=-RXONE/DCMPLX( S-DMZ**2 , DMAX1(DSIGN(DMZ*DWIDTH,S),RXZERO) )
      DZT=-RXONE/DCMPLX( T-DMZ**2 , DMAX1(DSIGN(DMZ*DWIDTH,T),RXZERO) )
C
      SV1= (DP2(0)+Q(0))*DV1(0) -(DP2(1)+Q(1))*DV1(1)
     &    -(DP2(2)+Q(2))*DV1(2) -(DP2(3)+Q(3))*DV1(3)
      SV2=-(DP1(0)+Q(0))*DV2(0) +(DP1(1)+Q(1))*DV2(1)
     &    +(DP1(2)+Q(2))*DV2(2) +(DP1(3)+Q(3))*DV2(3)
      SV3= (DP4(0)-Q(0))*DV3(0) -(DP4(1)-Q(1))*DV3(1)
     &    -(DP4(2)-Q(2))*DV3(2) -(DP4(3)-Q(3))*DV3(3)
      SV4=-(DP3(0)-Q(0))*DV4(0) +(DP3(1)-Q(1))*DV4(1)
     &    +(DP3(2)-Q(2))*DV4(2) +(DP3(3)-Q(3))*DV4(3)
C
      TV1= (DP4(0)+K(0))*DV1(0) -(DP4(1)+K(1))*DV1(1)
     &    -(DP4(2)+K(2))*DV1(2) -(DP4(3)+K(3))*DV1(3)
      TV2=-(DP3(0)-K(0))*DV2(0) +(DP3(1)-K(1))*DV2(1)
     &    +(DP3(2)-K(2))*DV2(2) +(DP3(3)-K(3))*DV2(3)
      TV3= (DP2(0)-K(0))*DV3(0) -(DP2(1)-K(1))*DV3(1)
     &    -(DP2(2)-K(2))*DV3(2) -(DP2(3)-K(3))*DV3(3)
      TV4=-(DP1(0)+K(0))*DV4(0) +(DP1(1)+K(1))*DV4(1)
     &    +(DP1(2)+K(2))*DV4(2) +(DP1(3)+K(3))*DV4(3)
C
      J12(0)=(DP1(0)-DP2(0))*V12 +SV1*DV2(0) +SV2*DV1(0)
      J12(1)=(DP1(1)-DP2(1))*V12 +SV1*DV2(1) +SV2*DV1(1)
      J12(2)=(DP1(2)-DP2(2))*V12 +SV1*DV2(2) +SV2*DV1(2)
      J12(3)=(DP1(3)-DP2(3))*V12 +SV1*DV2(3) +SV2*DV1(3)
      J34(0)=(DP3(0)-DP4(0))*V34 +SV3*DV4(0) +SV4*DV3(0)
      J34(1)=(DP3(1)-DP4(1))*V34 +SV3*DV4(1) +SV4*DV3(1)
      J34(2)=(DP3(2)-DP4(2))*V34 +SV3*DV4(2) +SV4*DV3(2)
      J34(3)=(DP3(3)-DP4(3))*V34 +SV3*DV4(3) +SV4*DV3(3)
C
      J14(0)=(DP1(0)-DP4(0))*V14 +TV1*DV4(0) +TV4*DV1(0)
      J14(1)=(DP1(1)-DP4(1))*V14 +TV1*DV4(1) +TV4*DV1(1)
      J14(2)=(DP1(2)-DP4(2))*V14 +TV1*DV4(2) +TV4*DV1(2)
      J14(3)=(DP1(3)-DP4(3))*V14 +TV1*DV4(3) +TV4*DV1(3)
      J32(0)=(DP3(0)-DP2(0))*V23 +TV3*DV2(0) +TV2*DV3(0)
      J32(1)=(DP3(1)-DP2(1))*V23 +TV3*DV2(1) +TV2*DV3(1)
      J32(2)=(DP3(2)-DP2(2))*V23 +TV3*DV2(2) +TV2*DV3(2)
      J32(3)=(DP3(3)-DP2(3))*V23 +TV3*DV2(3) +TV2*DV3(3)
C
      JS12=Q(0)*J12(0)-Q(1)*J12(1)-Q(2)*J12(2)-Q(3)*J12(3)
      JS34=Q(0)*J34(0)-Q(1)*J34(1)-Q(2)*J34(2)-Q(3)*J34(3)
      JS14=K(0)*J14(0)-K(1)*J14(1)-K(2)*J14(2)-K(3)*J14(3)
      JS32=K(0)*J32(0)-K(1)*J32(1)-K(2)*J32(2)-K(3)*J32(3)
C
      JS=J12(0)*J34(0)-J12(1)*J34(1)-J12(2)*J34(2)-J12(3)*J34(3)
      JT=J14(0)*J32(0)-J14(1)*J32(1)-J14(2)*J32(2)-J14(3)*J32(3)
C
      DVERTX = (V12*V34 +V14*V23 -RXTWO*V13*V24)*DGW2

C     &        +(dzs*dgwwz2+das*dgwwa2)*js -dzs*dgwwz2*js12*js34/dmz**2
C     &        +(dzt*dgwwz2+dat*dgwwa2)*jt -dzt*dgwwz2*js14*js32/dmz**2
C
      VERTEX = -DCMPLX( DVERTX )
C
      RETURN
      END
