* $Id$
    
C*********************************************************************  
    
      SUBROUTINE PYRESD_HIJING 
    
C...Allows resonances to decay (including parton showers for hadronic   
C...channels).  
      IMPLICIT DOUBLE PRECISION(D)  
#include "lujets_hijing.inc"
#include "ludat1_hijing.inc"
#include "ludat2_hijing.inc"
#include "ludat3_hijing.inc"
#include "pysubs_hijing.inc"
#include "pypars_hijing.inc"
#include "pyint1_hijing.inc"
#include "pyint2_hijing.inc"
#include "pyint4_hijing.inc"
      DIMENSION IREF(10,6),KDCY(2),KFL1(2),KFL2(2),NSD(2),ILIN(6),  
     &COUP(6,4),PK(6,4),PKK(6,6),CTHE(2),PHI(2),WDTP(0:40), 
     &WDTE(0:40,0:5)    
      COMPLEX FGK,HA(6,6),HC(6,6)   
    
C...The F, Xi and Xj functions of Gunion and Kunszt 
C...(Phys. Rev. D33, 665, plus errata from the authors).    
      FGK(I1,I2,I3,I4,I5,I6)=4.*HA(I1,I3)*HC(I2,I6)*(HA(I1,I5)* 
     &HC(I1,I4)+HA(I3,I5)*HC(I3,I4))    
      DIGK(DT,DU)=-4.*D34*D56+DT*(3.*DT+4.*DU)+DT**2*(DT*DU/(D34*D56)-  
     &2.*(1./D34+1./D56)*(DT+DU)+2.*(D34/D56+D56/D34))  
      DJGK(DT,DU)=8.*(D34+D56)**2-8.*(D34+D56)*(DT+DU)-6.*DT*DU-    
     &2.*DT*DU*(DT*DU/(D34*D56)-2.*(1./D34+1./D56)*(DT+DU)+ 
     &2.*(D34/D56+D56/D34)) 
    
C...Define initial two objects, initialize loop.    
      ISUB=MINT(1)  
      SH=VINT(44)   
      IREF(1,5)=0   
      IREF(1,6)=0   
      IF(ISET(ISUB).EQ.1.OR.ISET(ISUB).EQ.3) THEN   
        IREF(1,1)=MINT(84)+2+ISET(ISUB) 
        IREF(1,2)=0 
        IREF(1,3)=MINT(83)+6+ISET(ISUB) 
        IREF(1,4)=0 
      ELSEIF(ISET(ISUB).EQ.2.OR.ISET(ISUB).EQ.4) THEN   
        IREF(1,1)=MINT(84)+1+ISET(ISUB) 
        IREF(1,2)=MINT(84)+2+ISET(ISUB) 
        IREF(1,3)=MINT(83)+5+ISET(ISUB) 
        IREF(1,4)=MINT(83)+6+ISET(ISUB) 
      ENDIF 
      NP=1  
      IP=0  
  100 IP=IP+1   
      NINH=0    
    
C...Loop over one/two resonances; reset decay rates.    
      JTMAX=2   
      IF(IP.EQ.1.AND.(ISET(ISUB).EQ.1.OR.ISET(ISUB).EQ.3)) JTMAX=1  
      DO 140 JT=1,JTMAX 
      KDCY(JT)=0    
      KFL1(JT)=0    
      KFL2(JT)=0    
      NSD(JT)=IREF(IP,JT)   
      ID=IREF(IP,JT)    
      IF(ID.EQ.0) GOTO 140  
      KFA=IABS(K(ID,2)) 
      IF(KFA.LT.23.OR.KFA.GT.40) GOTO 140   
      IF(MDCY(KFA,1).NE.0) THEN 
        IF(ISUB.EQ.1.OR.ISUB.EQ.141) MINT(61)=1 
        CALL PYWIDT_HIJING(KFA,P(ID,5),WDTP,WDTE)  
        IF(KCHG(KFA,3).EQ.0) THEN   
          IPM=2 
        ELSE    
          IPM=(5+ISIGN(1,K(ID,2)))/2    
        ENDIF   
        IF(JTMAX.EQ.1.OR.IABS(K(IREF(IP,1),2)).NE.IABS(K(IREF(IP,2),2)))    
     &  THEN    
          I12=4 
        ELSE    
          IF(JT.EQ.1) I12=INT(4.5+RLU_HIJING(0))   
          I12=9-I12 
        ENDIF   
        RKFL=(WDTE(0,1)+WDTE(0,IPM)+WDTE(0,I12))*RLU_HIJING(0) 
        DO 120 I=1,MDCY(KFA,3)  
        IDC=I+MDCY(KFA,2)-1 
        KFL1(JT)=KFDP(IDC,1)*ISIGN(1,K(ID,2))   
        KFL2(JT)=KFDP(IDC,2)*ISIGN(1,K(ID,2))   
        RKFL=RKFL-(WDTE(I,1)+WDTE(I,IPM)+WDTE(I,I12))   
        IF(RKFL.LE.0.) GOTO 130 
  120   CONTINUE    
  130   CONTINUE    
      ENDIF 
    
C...Summarize result on decay channel chosen.   
      IF((KFA.EQ.23.OR.KFA.EQ.24).AND.KFL1(JT).EQ.0) NINH=NINH+1    
      IF(KFL1(JT).EQ.0) GOTO 140    
      KDCY(JT)=2    
      IF(IABS(KFL1(JT)).LE.10.OR.KFL1(JT).EQ.21) KDCY(JT)=1 
      IF((IABS(KFL1(JT)).GE.23.AND.IABS(KFL1(JT)).LE.25).OR.    
     &(IABS(KFL1(JT)).EQ.37)) KDCY(JT)=3    
      NSD(JT)=N 
    
C...Fill decay products, prepared for parton showers for quarks.    
      IF(KDCY(JT).EQ.1) THEN    
        CALL LU2ENT_HIJING(-(N+1),KFL1(JT),KFL2(JT),P(ID,5))   
      ELSE  
        CALL LU2ENT_HIJING(N+1,KFL1(JT),KFL2(JT),P(ID,5))  
      ENDIF 
      IF(JTMAX.EQ.1) THEN   
         CTHE(JT)=VINT(13)+(VINT(33)-VINT(13)+VINT(34)-VINT(14))
     $        *RLU_HIJING(0)  
        IF(CTHE(JT).GT.VINT(33)) CTHE(JT)=CTHE(JT)+VINT(14)-VINT(33)    
        PHI(JT)=VINT(24)    
      ELSE  
        CTHE(JT)=2.*RLU_HIJING(0)-1.   
        PHI(JT)=PARU(2)*RLU_HIJING(0)  
      ENDIF 
  140 CONTINUE  
      IF(MINT(3).EQ.1.AND.IP.EQ.1) THEN 
        MINT(25)=KFL1(1)    
        MINT(26)=KFL2(1)    
      ENDIF 
      IF(JTMAX.EQ.1.AND.KDCY(1).EQ.0) GOTO 530  
      IF(JTMAX.EQ.2.AND.KDCY(1).EQ.0.AND.KDCY(2).EQ.0) GOTO 530 
      IF(MSTP(45).LE.0.OR.IREF(IP,2).EQ.0.OR.NINH.GE.1) GOTO 500    
      IF(K(IREF(1,1),2).EQ.25.AND.IP.EQ.1) GOTO 500 
      IF(K(IREF(1,1),2).EQ.25.AND.KDCY(1)*KDCY(2).EQ.0) GOTO 500    
    
C...Order incoming partons and outgoing resonances. 
      ILIN(1)=MINT(84)+1    
      IF(K(MINT(84)+1,2).GT.0) ILIN(1)=MINT(84)+2   
      IF(K(ILIN(1),2).EQ.21) ILIN(1)=2*MINT(84)+3-ILIN(1)   
      ILIN(2)=2*MINT(84)+3-ILIN(1)  
      IMIN=1    
      IF(IREF(IP,5).EQ.25) IMIN=3   
      IMAX=2    
      IORD=1    
      IF(K(IREF(IP,1),2).EQ.23) IORD=2  
      IF(K(IREF(IP,1),2).EQ.24.AND.K(IREF(IP,2),2).EQ.-24) IORD=2   
      IF(IABS(K(IREF(IP,IORD),2)).EQ.25) IORD=3-IORD    
      IF(KDCY(IORD).EQ.0) IORD=3-IORD   
    
C...Order decay products of resonances. 
      DO 390 JT=IORD,3-IORD,3-2*IORD    
      IF(KDCY(JT).EQ.0) THEN    
        ILIN(IMAX+1)=NSD(JT)    
        IMAX=IMAX+1 
      ELSEIF(K(NSD(JT)+1,2).GT.0) THEN  
        ILIN(IMAX+1)=N+2*JT-1   
        ILIN(IMAX+2)=N+2*JT 
        IMAX=IMAX+2 
        K(N+2*JT-1,2)=K(NSD(JT)+1,2)    
        K(N+2*JT,2)=K(NSD(JT)+2,2)  
      ELSE  
        ILIN(IMAX+1)=N+2*JT 
        ILIN(IMAX+2)=N+2*JT-1   
        IMAX=IMAX+2 
        K(N+2*JT-1,2)=K(NSD(JT)+1,2)    
        K(N+2*JT,2)=K(NSD(JT)+2,2)  
      ENDIF 
  390 CONTINUE  
    
C...Find charge, isospin, left- and righthanded couplings.  
      XW=PARU(102)  
      DO 410 I=IMIN,IMAX    
      DO 400 J=1,4  
  400 COUP(I,J)=0.  
      KFA=IABS(K(ILIN(I),2))    
      IF(KFA.GT.20) GOTO 410    
      COUP(I,1)=LUCHGE_HIJING(KFA)/3.  
      COUP(I,2)=(-1)**MOD(KFA,2)    
      COUP(I,4)=-2.*COUP(I,1)*XW    
      COUP(I,3)=COUP(I,2)+COUP(I,4) 
  410 CONTINUE  
      SQMZ=PMAS(23,1)**2    
      GZMZ=PMAS(23,1)*PMAS(23,2)    
      SQMW=PMAS(24,1)**2    
      GZMW=PMAS(24,1)*PMAS(24,2)    
      SQMZP=PMAS(32,1)**2   
      GZMZP=PMAS(32,1)*PMAS(32,2)   
    
C...Select random angles; construct massless four-vectors.  
  420 DO 430 I=N+1,N+4  
      K(I,1)=1  
      DO 430 J=1,5  
  430 P(I,J)=0. 
      DO 440 JT=1,JTMAX 
      IF(KDCY(JT).EQ.0) GOTO 440    
      ID=IREF(IP,JT)    
      P(N+2*JT-1,3)=0.5*P(ID,5) 
      P(N+2*JT-1,4)=0.5*P(ID,5) 
      P(N+2*JT,3)=-0.5*P(ID,5)  
      P(N+2*JT,4)=0.5*P(ID,5)   
      CTHE(JT)=2.*RLU_HIJING(0)-1. 
      PHI(JT)=PARU(2)*RLU_HIJING(0)    
      CALL LUDBRB_HIJING(N+2*JT-1,N+2*JT,ACOS(CTHE(JT)),PHI(JT),   
     &DBLE(P(ID,1)/P(ID,4)),DBLE(P(ID,2)/P(ID,4)),DBLE(P(ID,3)/P(ID,4)))    
  440 CONTINUE  
    
C...Store incoming and outgoing momenta, with random rotation to    
C...avoid accidental zeroes in HA expressions.  
      DO 450 I=1,IMAX   
      K(N+4+I,1)=1  
      P(N+4+I,4)=SQRT(P(ILIN(I),1)**2+P(ILIN(I),2)**2+P(ILIN(I),3)**2+  
     &P(ILIN(I),5)**2)  
      P(N+4+I,5)=P(ILIN(I),5)   
      DO 450 J=1,3  
  450 P(N+4+I,J)=P(ILIN(I),J)   
      THERR=ACOS(2.*RLU_HIJING(0)-1.)  
      PHIRR=PARU(2)*RLU_HIJING(0)  
      CALL LUDBRB_HIJING(N+5,N+4+IMAX,THERR,PHIRR,0D0,0D0,0D0) 
      DO 460 I=1,IMAX   
      DO 460 J=1,4  
  460 PK(I,J)=P(N+4+I,J)    
    
C...Calculate internal products.    
      IF(ISUB.EQ.22.OR.ISUB.EQ.23.OR.ISUB.EQ.25) THEN   
        DO 470 I1=IMIN,IMAX-1   
        DO 470 I2=I1+1,IMAX 
        HA(I1,I2)=SQRT((PK(I1,4)-PK(I1,3))*(PK(I2,4)+PK(I2,3))/ 
     &  (1E-20+PK(I1,1)**2+PK(I1,2)**2))*CMPLX(PK(I1,1),PK(I1,2))-  
     &  SQRT((PK(I1,4)+PK(I1,3))*(PK(I2,4)-PK(I2,3))/   
     &  (1E-20+PK(I2,1)**2+PK(I2,2)**2))*CMPLX(PK(I2,1),PK(I2,2))   
        HC(I1,I2)=CONJG(HA(I1,I2))  
        IF(I1.LE.2) HA(I1,I2)=CMPLX(0.,1.)*HA(I1,I2)    
        IF(I1.LE.2) HC(I1,I2)=CMPLX(0.,1.)*HC(I1,I2)    
        HA(I2,I1)=-HA(I1,I2)    
  470   HC(I2,I1)=-HC(I1,I2)    
      ENDIF 
      DO 480 I=1,2  
      DO 480 J=1,4  
  480 PK(I,J)=-PK(I,J)  
      DO 490 I1=IMIN,IMAX-1 
      DO 490 I2=I1+1,IMAX   
      PKK(I1,I2)=2.*(PK(I1,4)*PK(I2,4)-PK(I1,1)*PK(I2,1)-   
     &PK(I1,2)*PK(I2,2)-PK(I1,3)*PK(I2,3))  
  490 PKK(I2,I1)=PKK(I1,I2) 
    
      IF(IREF(IP,5).EQ.25) THEN 
C...Angular weight for H0 -> Z0 + Z0 or W+ + W- -> 4 quarks/leptons 
        WT=16.*PKK(3,5)*PKK(4,6)    
        IF(IP.EQ.1) WTMAX=SH**2 
        IF(IP.GE.2) WTMAX=P(IREF(IP,6),5)**4    
    
      ELSEIF(ISUB.EQ.1) THEN    
        IF(KFA.NE.37) THEN  
C...Angular weight for gamma*/Z0 -> 2 quarks/leptons    
          EI=KCHG(IABS(MINT(15)),1)/3.  
          AI=SIGN(1.,EI+0.1)    
          VI=AI-4.*EI*XW    
          EF=KCHG(KFA,1)/3. 
          AF=SIGN(1.,EF+0.1)    
          VF=AF-4.*EF*XW    
          GG=1. 
          GZ=1./(8.*XW*(1.-XW))*SH*(SH-SQMZ)/((SH-SQMZ)**2+GZMZ**2) 
          ZZ=1./(16.*XW*(1.-XW))**2*SH**2/((SH-SQMZ)**2+GZMZ**2)    
          IF(MSTP(43).EQ.1) THEN    
C...Only gamma* production included 
            GZ=0.   
            ZZ=0.   
          ELSEIF(MSTP(43).EQ.2) THEN    
C...Only Z0 production included 
            GG=0.   
            GZ=0.   
          ENDIF 
          ASYM=2.*(EI*AI*GZ*EF*AF+4.*VI*AI*ZZ*VF*AF)/(EI**2*GG*EF**2+   
     &    EI*VI*GZ*EF*VF+(VI**2+AI**2)*ZZ*(VF**2+AF**2))    
          WT=1.+ASYM*CTHE(JT)+CTHE(JT)**2   
          WTMAX=2.+ABS(ASYM)    
        ELSE    
C...Angular weight for gamma*/Z0 -> H+ + H- 
          WT=1.-CTHE(JT)**2 
          WTMAX=1.  
        ENDIF   
    
      ELSEIF(ISUB.EQ.2) THEN    
C...Angular weight for W+/- -> 2 quarks/leptons 
        WT=(1.+CTHE(JT))**2 
        WTMAX=4.    
    
      ELSEIF(ISUB.EQ.15.OR.ISUB.EQ.19) THEN 
C...Angular weight for f + fb -> gluon/gamma + Z0 ->    
C...-> gluon/gamma + 2 quarks/leptons   
        WT=((COUP(1,3)*COUP(3,3))**2+(COUP(1,4)*COUP(3,4))**2)* 
     &  (PKK(1,3)**2+PKK(2,4)**2)+((COUP(1,3)*COUP(3,4))**2+    
     &  (COUP(1,4)*COUP(3,3))**2)*(PKK(1,4)**2+PKK(2,3)**2) 
        WTMAX=(COUP(1,3)**2+COUP(1,4)**2)*(COUP(3,3)**2+COUP(3,4)**2)*  
     &  ((PKK(1,3)+PKK(1,4))**2+(PKK(2,3)+PKK(2,4))**2) 
    
      ELSEIF(ISUB.EQ.16.OR.ISUB.EQ.20) THEN 
C...Angular weight for f + fb' -> gluon/gamma + W+/- -> 
C...-> gluon/gamma + 2 quarks/leptons   
        WT=PKK(1,3)**2+PKK(2,4)**2  
        WTMAX=(PKK(1,3)+PKK(1,4))**2+(PKK(2,3)+PKK(2,4))**2 
    
      ELSEIF(ISUB.EQ.22) THEN   
C...Angular weight for f + fb -> Z0 + Z0 -> 4 quarks/leptons    
        S34=P(IREF(IP,IORD),5)**2   
        S56=P(IREF(IP,3-IORD),5)**2 
        TI=PKK(1,3)+PKK(1,4)+S34    
        UI=PKK(1,5)+PKK(1,6)+S56    
        WT=COUP(1,3)**4*((COUP(3,3)*COUP(5,3)*ABS(FGK(1,2,3,4,5,6)/ 
     &  TI+FGK(1,2,5,6,3,4)/UI))**2+(COUP(3,4)*COUP(5,3)*ABS(   
     &  FGK(1,2,4,3,5,6)/TI+FGK(1,2,5,6,4,3)/UI))**2+(COUP(3,3)*    
     &  COUP(5,4)*ABS(FGK(1,2,3,4,6,5)/TI+FGK(1,2,6,5,3,4)/UI))**2+ 
     &  (COUP(3,4)*COUP(5,4)*ABS(FGK(1,2,4,3,6,5)/TI+FGK(1,2,6,5,4,3)/  
     &  UI))**2)+COUP(1,4)**4*((COUP(3,3)*COUP(5,3)*ABS(    
     &  FGK(2,1,5,6,3,4)/TI+FGK(2,1,3,4,5,6)/UI))**2+(COUP(3,4)*    
     &  COUP(5,3)*ABS(FGK(2,1,6,5,3,4)/TI+FGK(2,1,3,4,6,5)/UI))**2+ 
     &  (COUP(3,3)*COUP(5,4)*ABS(FGK(2,1,5,6,4,3)/TI+FGK(2,1,4,3,5,6)/  
     &  UI))**2+(COUP(3,4)*COUP(5,4)*ABS(FGK(2,1,6,5,4,3)/TI+   
     &  FGK(2,1,4,3,6,5)/UI))**2)   
        WTMAX=4.*S34*S56*(COUP(1,3)**4+COUP(1,4)**4)*(COUP(3,3)**2+ 
     &  COUP(3,4)**2)*(COUP(5,3)**2+COUP(5,4)**2)*4.*(TI/UI+UI/TI+  
     &  2.*SH*(S34+S56)/(TI*UI)-S34*S56*(1./TI**2+1./UI**2))    
    
      ELSEIF(ISUB.EQ.23) THEN   
C...Angular weight for f + fb' -> Z0 + W +/- -> 4 quarks/leptons    
        D34=P(IREF(IP,IORD),5)**2   
        D56=P(IREF(IP,3-IORD),5)**2 
        DT=PKK(1,3)+PKK(1,4)+D34    
        DU=PKK(1,5)+PKK(1,6)+D56    
        CAWZ=COUP(2,3)/SNGL(DT)-2.*(1.-XW)*COUP(1,2)/(SH-SQMW)  
        CBWZ=COUP(1,3)/SNGL(DU)+2.*(1.-XW)*COUP(1,2)/(SH-SQMW)  
        WT=COUP(5,3)**2*ABS(CAWZ*FGK(1,2,3,4,5,6)+CBWZ* 
     &  FGK(1,2,5,6,3,4))**2+COUP(5,4)**2*ABS(CAWZ* 
     &  FGK(1,2,3,4,6,5)+CBWZ*FGK(1,2,6,5,3,4))**2  
        WTMAX=4.*D34*D56*(COUP(5,3)**2+COUP(5,4)**2)*(CAWZ**2*  
     &  DIGK(DT,DU)+CBWZ**2*DIGK(DU,DT)+CAWZ*CBWZ*DJGK(DT,DU))  
    
      ELSEIF(ISUB.EQ.24) THEN   
C...Angular weight for f + fb -> Z0 + H0 -> 2 quarks/leptons + H0   
        WT=((COUP(1,3)*COUP(3,3))**2+(COUP(1,4)*COUP(3,4))**2)* 
     &  PKK(1,3)*PKK(2,4)+((COUP(1,3)*COUP(3,4))**2+(COUP(1,4)* 
     &  COUP(3,3))**2)*PKK(1,4)*PKK(2,3)    
        WTMAX=(COUP(1,3)**2+COUP(1,4)**2)*(COUP(3,3)**2+COUP(3,4)**2)*  
     &  (PKK(1,3)+PKK(1,4))*(PKK(2,3)+PKK(2,4)) 
    
      ELSEIF(ISUB.EQ.25) THEN   
C...Angular weight for f + fb -> W+ + W- -> 4 quarks/leptons    
        D34=P(IREF(IP,IORD),5)**2   
        D56=P(IREF(IP,3-IORD),5)**2 
        DT=PKK(1,3)+PKK(1,4)+D34    
        DU=PKK(1,5)+PKK(1,6)+D56    
        CDWW=(COUP(1,3)*SQMZ/(SH-SQMZ)+COUP(1,2))/SH    
        CAWW=CDWW+0.5*(COUP(1,2)+1.)/SNGL(DT)   
        CBWW=CDWW+0.5*(COUP(1,2)-1.)/SNGL(DU)   
        CCWW=COUP(1,4)*SQMZ/(SH-SQMZ)/SH    
        WT=ABS(CAWW*FGK(1,2,3,4,5,6)-CBWW*FGK(1,2,5,6,3,4))**2+ 
     &  CCWW**2*ABS(FGK(2,1,5,6,3,4)-FGK(2,1,3,4,5,6))**2   
        WTMAX=4.*D34*D56*(CAWW**2*DIGK(DT,DU)+CBWW**2*DIGK(DU,DT)-CAWW* 
     &  CBWW*DJGK(DT,DU)+CCWW**2*(DIGK(DT,DU)+DIGK(DU,DT)-DJGK(DT,DU))) 
    
      ELSEIF(ISUB.EQ.26) THEN   
C...Angular weight for f + fb' -> W+/- + H0 -> 2 quarks/leptons + H0    
        WT=PKK(1,3)*PKK(2,4)    
        WTMAX=(PKK(1,3)+PKK(1,4))*(PKK(2,3)+PKK(2,4))   
    
      ELSEIF(ISUB.EQ.30) THEN   
C...Angular weight for f + g -> f + Z0 -> f + 2 quarks/leptons  
        IF(K(ILIN(1),2).GT.0) WT=((COUP(1,3)*COUP(3,3))**2+ 
     &  (COUP(1,4)*COUP(3,4))**2)*(PKK(1,4)**2+PKK(3,5)**2)+    
     &  ((COUP(1,3)*COUP(3,4))**2+(COUP(1,4)*COUP(3,3))**2)*    
     &  (PKK(1,3)**2+PKK(4,5)**2)   
        IF(K(ILIN(1),2).LT.0) WT=((COUP(1,3)*COUP(3,3))**2+ 
     &  (COUP(1,4)*COUP(3,4))**2)*(PKK(1,3)**2+PKK(4,5)**2)+    
     &  ((COUP(1,3)*COUP(3,4))**2+(COUP(1,4)*COUP(3,3))**2)*    
     &  (PKK(1,4)**2+PKK(3,5)**2)   
        WTMAX=(COUP(1,3)**2+COUP(1,4)**2)*(COUP(3,3)**2+COUP(3,4)**2)*  
     &  ((PKK(1,3)+PKK(1,4))**2+(PKK(3,5)+PKK(4,5))**2) 
    
      ELSEIF(ISUB.EQ.31) THEN   
C...Angular weight for f + g -> f' + W+/- -> f' + 2 quarks/leptons  
        IF(K(ILIN(1),2).GT.0) WT=PKK(1,4)**2+PKK(3,5)**2    
        IF(K(ILIN(1),2).LT.0) WT=PKK(1,3)**2+PKK(4,5)**2    
        WTMAX=(PKK(1,3)+PKK(1,4))**2+(PKK(3,5)+PKK(4,5))**2 
    
      ELSEIF(ISUB.EQ.141) THEN  
C...Angular weight for gamma*/Z0/Z'0 -> 2 quarks/leptons    
        EI=KCHG(IABS(MINT(15)),1)/3.    
        AI=SIGN(1.,EI+0.1)  
        VI=AI-4.*EI*XW  
        API=SIGN(1.,EI+0.1) 
        VPI=API-4.*EI*XW    
        EF=KCHG(KFA,1)/3.   
        AF=SIGN(1.,EF+0.1)  
        VF=AF-4.*EF*XW  
        APF=SIGN(1.,EF+0.1) 
        VPF=APF-4.*EF*XW    
        GG=1.   
        GZ=1./(8.*XW*(1.-XW))*SH*(SH-SQMZ)/((SH-SQMZ)**2+GZMZ**2)   
        GZP=1./(8.*XW*(1.-XW))*SH*(SH-SQMZP)/((SH-SQMZP)**2+GZMZP**2)   
        ZZ=1./(16.*XW*(1.-XW))**2*SH**2/((SH-SQMZ)**2+GZMZ**2)  
        ZZP=2./(16.*XW*(1.-XW))**2* 
     &  SH**2*((SH-SQMZ)*(SH-SQMZP)+GZMZ*GZMZP)/    
     &  (((SH-SQMZ)**2+GZMZ**2)*((SH-SQMZP)**2+GZMZP**2))   
        ZPZP=1./(16.*XW*(1.-XW))**2*SH**2/((SH-SQMZP)**2+GZMZP**2)  
        IF(MSTP(44).EQ.1) THEN  
C...Only gamma* production included 
          GZ=0. 
          GZP=0.    
          ZZ=0. 
          ZZP=0.    
          ZPZP=0.   
        ELSEIF(MSTP(44).EQ.2) THEN  
C...Only Z0 production included 
          GG=0. 
          GZ=0. 
          GZP=0.    
          ZZP=0.    
          ZPZP=0.   
        ELSEIF(MSTP(44).EQ.3) THEN  
C...Only Z'0 production included    
          GG=0. 
          GZ=0. 
          GZP=0.    
          ZZ=0. 
          ZZP=0.    
        ELSEIF(MSTP(44).EQ.4) THEN  
C...Only gamma*/Z0 production included  
          GZP=0.    
          ZZP=0.    
          ZPZP=0.   
        ELSEIF(MSTP(44).EQ.5) THEN  
C...Only gamma*/Z'0 production included 
          GZ=0. 
          ZZ=0. 
          ZZP=0.    
        ELSEIF(MSTP(44).EQ.6) THEN  
C...Only Z0/Z'0 production included 
          GG=0. 
          GZ=0. 
          GZP=0.    
        ENDIF   
        ASYM=2.*(EI*AI*GZ*EF*AF+EI*API*GZP*EF*APF+4.*VI*AI*ZZ*VF*AF+    
     &  (VI*API+VPI*AI)*ZZP*(VF*APF+VPF*AF)+4.*VPI*API*ZPZP*VPF*APF)/   
     &  (EI**2*GG*EF**2+EI*VI*GZ*EF*VF+EI*VPI*GZP*EF*VPF+   
     &  (VI**2+AI**2)*ZZ*(VF**2+AF**2)+(VI*VPI+AI*API)*ZZP* 
     &  (VF*VPF+AF*APF)+(VPI**2+API**2)*ZPZP*(VPF**2+APF**2))   
        WT=1.+ASYM*CTHE(JT)+CTHE(JT)**2 
        WTMAX=2.+ABS(ASYM)  
    
      ELSE  
        WT=1.   
        WTMAX=1.    
      ENDIF 
C...Obtain correct angular distribution by rejection techniques.    
      IF(WT.LT.RLU_HIJING(0)*WTMAX) GOTO 420   
    
C...Construct massive four-vectors using angles chosen. Mark decayed    
C...resonances, add documentation lines. Shower evolution.  
  500 DO 520 JT=1,JTMAX 
      IF(KDCY(JT).EQ.0) GOTO 520    
      ID=IREF(IP,JT)    
      CALL LUDBRB_HIJING(NSD(JT)+1,NSD(JT)+2,ACOS(CTHE(JT)),PHI(JT),   
     &DBLE(P(ID,1)/P(ID,4)),DBLE(P(ID,2)/P(ID,4)),DBLE(P(ID,3)/P(ID,4)))    
      K(ID,1)=K(ID,1)+10    
      K(ID,4)=NSD(JT)+1 
      K(ID,5)=NSD(JT)+2 
      IDOC=MINT(83)+MINT(4) 
      DO 510 I=NSD(JT)+1,NSD(JT)+2  
      MINT(4)=MINT(4)+1 
      I1=MINT(83)+MINT(4)   
      K(I,3)=I1 
      K(I1,1)=21    
      K(I1,2)=K(I,2)    
      K(I1,3)=IREF(IP,JT+2) 
      DO 510 J=1,5  
  510 P(I1,J)=P(I,J)    
      IF(JTMAX.EQ.1) THEN   
        MINT(7)=MINT(83)+6+2*ISET(ISUB) 
        MINT(8)=MINT(83)+7+2*ISET(ISUB) 
      ENDIF 
      IF(MSTP(71).GE.1.AND.KDCY(JT).EQ.1) CALL LUSHOW_HIJING(NSD(JT)+1,    
     &NSD(JT)+2,P(ID,5))    
    
C...Check if new resonances were produced, loop back if needed. 
      IF(KDCY(JT).NE.3) GOTO 520    
      NP=NP+1   
      IREF(NP,1)=NSD(JT)+1  
      IREF(NP,2)=NSD(JT)+2  
      IREF(NP,3)=IDOC+1 
      IREF(NP,4)=IDOC+2 
      IREF(NP,5)=K(IREF(IP,JT),2)   
      IREF(NP,6)=IREF(IP,JT)    
  520 CONTINUE  
  530 IF(IP.LT.NP) GOTO 100 
    
      RETURN    
      END   
