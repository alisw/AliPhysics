* $Id$
    
C***********************************************************************    
    
      SUBROUTINE PYSIGH_HIJING(NCHN,SIGS)  
    
C...Differential matrix elements for all included subprocesses. 
C...Note that what is coded is (disregarding the COMFAC factor) 
C...1) for 2 -> 1 processes: s-hat/pi*d(sigma-hat), where,  
C...when d(sigma-hat) is given in the zero-width limit, the delta   
C...function in tau is replaced by a Breit-Wigner:  
C...1/pi*(s*m_res*Gamma_res)/((s*tau-m_res^2)^2+(m_res*Gamma_res)^2);   
C...2) for 2 -> 2 processes: (s-hat)**2/pi*d(sigma-hat)/d(t-hat);   
C...i.e., dimensionless quantities. COMFAC contains the factor  
C...pi/s and the conversion factor from GeV^-2 to mb.   
#include "ludat1_hijing.inc"
#include "ludat2_hijing.inc"
#include "ludat3_hijing.inc"
#include "pysubs_hijing.inc"
#include "pypars_hijing.inc"
#include "pyint1_hijing.inc"
#include "pyint2_hijing.inc"
#include "pyint3_hijing.inc"
#include "pyint4_hijing.inc"
#include "pyint5_hijing.inc"
      DIMENSION X(2),XPQ(-6:6),KFAC(2,-40:40),WDTP(0:40),WDTE(0:40,0:5) 
    
C...Reset number of channels and cross-section. 
      NCHN=0    
      SIGS=0.   
    
C...Read kinematical variables and limits.  
      ISUB=MINT(1)  
      TAUMIN=VINT(11)   
      YSTMIN=VINT(12)   
      CTNMIN=VINT(13)   
      CTPMIN=VINT(14)   
      XT2MIN=VINT(15)   
      TAUPMN=VINT(16)   
      TAU=VINT(21)  
      YST=VINT(22)  
      CTH=VINT(23)  
      XT2=VINT(25)  
      TAUP=VINT(26) 
      TAUMAX=VINT(31)   
      YSTMAX=VINT(32)   
      CTNMAX=VINT(33)   
      CTPMAX=VINT(34)   
      XT2MAX=VINT(35)   
      TAUPMX=VINT(36)   
    
C...Derive kinematical quantities.  
      IF(ISET(ISUB).LE.2.OR.ISET(ISUB).EQ.5) THEN   
        X(1)=SQRT(TAU)*EXP(YST) 
        X(2)=SQRT(TAU)*EXP(-YST)    
      ELSE  
        X(1)=SQRT(TAUP)*EXP(YST)    
        X(2)=SQRT(TAUP)*EXP(-YST)   
      ENDIF 
      IF(MINT(43).EQ.4.AND.ISET(ISUB).GE.1.AND. 
     &(X(1).GT.0.999.OR.X(2).GT.0.999)) RETURN  
      SH=TAU*VINT(2)    
      SQM3=VINT(63) 
      SQM4=VINT(64) 
      RM3=SQM3/SH   
      RM4=SQM4/SH   
      BE34=SQRT((1.-RM3-RM4)**2-4.*RM3*RM4) 
      RPTS=4.*VINT(71)**2/SH    
      BE34L=SQRT(MAX(0.,(1.-RM3-RM4)**2-4.*RM3*RM4-RPTS))   
      RM34=2.*RM3*RM4   
      RSQM=1.+RM34  
      RTHM=(4.*RM3*RM4+RPTS)/(1.-RM3-RM4+BE34L) 
      TH=-0.5*SH*MAX(RTHM,1.-RM3-RM4-BE34*CTH)  
      UH=-0.5*SH*MAX(RTHM,1.-RM3-RM4+BE34*CTH)  
      SQPTH=0.25*SH*BE34**2*(1.-CTH**2) 
      SH2=SH**2 
      TH2=TH**2 
      UH2=UH**2 
    
C...Choice of Q2 scale. 
      IF(ISET(ISUB).EQ.1.OR.ISET(ISUB).EQ.3) THEN   
        Q2=SH   
      ELSEIF(MOD(ISET(ISUB),2).EQ.0.OR.ISET(ISUB).EQ.5) THEN    
        IF(MSTP(32).EQ.1) THEN  
          Q2=2.*SH*TH*UH/(SH**2+TH**2+UH**2)    
        ELSEIF(MSTP(32).EQ.2) THEN  
          Q2=SQPTH+0.5*(SQM3+SQM4)  
        ELSEIF(MSTP(32).EQ.3) THEN  
          Q2=MIN(-TH,-UH)   
        ELSEIF(MSTP(32).EQ.4) THEN  
          Q2=SH 
        ENDIF   
        IF(ISET(ISUB).EQ.5.AND.MSTP(82).GE.2) Q2=Q2+PARP(82)**2 
      ENDIF 
    
C...Store derived kinematical quantities.   
      VINT(41)=X(1) 
      VINT(42)=X(2) 
      VINT(44)=SH   
      VINT(43)=SQRT(SH) 
      VINT(45)=TH   
      VINT(46)=UH   
      VINT(48)=SQPTH    
      VINT(47)=SQRT(SQPTH)  
      VINT(50)=TAUP*VINT(2) 
      VINT(49)=SQRT(MAX(0.,VINT(50)))   
      VINT(52)=Q2   
      VINT(51)=SQRT(Q2) 
    
C...Calculate parton structure functions.   
      IF(ISET(ISUB).LE.0) GOTO 145  
      IF(MINT(43).GE.2) THEN    
        Q2SF=Q2 
        IF(ISET(ISUB).EQ.3.OR.ISET(ISUB).EQ.4) THEN 
          Q2SF=PMAS(23,1)**2    
          IF(ISUB.EQ.8.OR.ISUB.EQ.76.OR.ISUB.EQ.77) Q2SF=PMAS(24,1)**2  
        ENDIF   
        DO 100 I=3-MINT(41),MINT(42)    
        XSF=X(I)    
        IF(ISET(ISUB).EQ.5) XSF=X(I)/VINT(142+I)    
        CALL PYSTFU_HIJING(MINT(10+I),XSF,Q2SF,XPQ,I)    
        DO 100 KFL=-6,6 
  100   XSFX(I,KFL)=XPQ(KFL)
      ENDIF 
    
C...Calculate alpha_strong and K-factor.    
      IF(MSTP(33).NE.3) AS=ULALPS_HIJING(Q2)   
      FACK=1.   
      FACA=1.   
      IF(MSTP(33).EQ.1) THEN    
        FACK=PARP(31)   
      ELSEIF(MSTP(33).EQ.2) THEN    
        FACK=PARP(31)   
        FACA=PARP(32)/PARP(31)  
      ELSEIF(MSTP(33).EQ.3) THEN    
        Q2AS=PARP(33)*Q2    
        IF(ISET(ISUB).EQ.5.AND.MSTP(82).GE.2) Q2AS=Q2AS+    
     &  PARU(112)*PARP(82)  
        AS=ULALPS_HIJING(Q2AS) 
      ENDIF 
      RADC=1.+AS/PARU(1)    
    
C...Set flags for allowed reacting partons/leptons. 
      DO 130 I=1,2  
      DO 110 J=-40,40   
  110 KFAC(I,J)=0   
      IF(MINT(40+I).EQ.1) THEN  
        KFAC(I,MINT(10+I))=1    
      ELSE  
        DO 120 J=-40,40 
        KFAC(I,J)=KFIN(I,J) 
        IF(ABS(J).GT.MSTP(54).AND.J.NE.21) KFAC(I,J)=0  
        IF(ABS(J).LE.6) THEN    
          IF(XSFX(I,J).LT.1.E-10) KFAC(I,J)=0   
        ELSEIF(J.EQ.21) THEN    
          IF(XSFX(I,0).LT.1.E-10) KFAC(I,21)=0  
        ENDIF   
  120   CONTINUE    
      ENDIF 
  130 CONTINUE  
    
C...Lower and upper limit for flavour loops.    
      MIN1=0    
      MAX1=0    
      MIN2=0    
      MAX2=0    
      DO 140 J=-20,20   
      IF(KFAC(1,-J).EQ.1) MIN1=-J   
      IF(KFAC(1,J).EQ.1) MAX1=J 
      IF(KFAC(2,-J).EQ.1) MIN2=-J   
      IF(KFAC(2,J).EQ.1) MAX2=J 
  140 CONTINUE  
      MINA=MIN(MIN1,MIN2)   
      MAXA=MAX(MAX1,MAX2)   
    
C...Common conversion factors (including Jacobian) for subprocesses.    
      SQMZ=PMAS(23,1)**2    
      GMMZ=PMAS(23,1)*PMAS(23,2)    
      SQMW=PMAS(24,1)**2    
      GMMW=PMAS(24,1)*PMAS(24,2)    
      SQMH=PMAS(25,1)**2    
      GMMH=PMAS(25,1)*PMAS(25,2)    
      SQMZP=PMAS(32,1)**2   
      GMMZP=PMAS(32,1)*PMAS(32,2)   
      SQMHC=PMAS(37,1)**2   
      GMMHC=PMAS(37,1)*PMAS(37,2)   
      SQMR=PMAS(40,1)**2    
      GMMR=PMAS(40,1)*PMAS(40,2)    
      AEM=PARU(101) 
      XW=PARU(102)  
    
C...Phase space integral in tau and y*. 
      COMFAC=PARU(1)*PARU(5)/VINT(2)    
      IF(MINT(43).EQ.4) COMFAC=COMFAC*FACK  
      IF((MINT(43).GE.2.OR.ISET(ISUB).EQ.3.OR.ISET(ISUB).EQ.4).AND. 
     &ISET(ISUB).NE.5) THEN 
        ATAU0=LOG(TAUMAX/TAUMIN)    
        ATAU1=(TAUMAX-TAUMIN)/(TAUMAX*TAUMIN)   
        H1=COEF(ISUB,1)+(ATAU0/ATAU1)*COEF(ISUB,2)/TAU  
        IF(MINT(72).GE.1) THEN  
          TAUR1=VINT(73)    
          GAMR1=VINT(74)    
          ATAU2=LOG(TAUMAX/TAUMIN*(TAUMIN+TAUR1)/(TAUMAX+TAUR1))/TAUR1  
          ATAU3=(ATAN((TAUMAX-TAUR1)/GAMR1)-ATAN((TAUMIN-TAUR1)/GAMR1))/    
     &    GAMR1 
          H1=H1+(ATAU0/ATAU2)*COEF(ISUB,3)/(TAU+TAUR1)+ 
     &    (ATAU0/ATAU3)*COEF(ISUB,4)*TAU/((TAU-TAUR1)**2+GAMR1**2)  
        ENDIF   
        IF(MINT(72).EQ.2) THEN  
          TAUR2=VINT(75)    
          GAMR2=VINT(76)    
          ATAU4=LOG(TAUMAX/TAUMIN*(TAUMIN+TAUR2)/(TAUMAX+TAUR2))/TAUR2  
          ATAU5=(ATAN((TAUMAX-TAUR2)/GAMR2)-ATAN((TAUMIN-TAUR2)/GAMR2))/    
     &    GAMR2 
          H1=H1+(ATAU0/ATAU4)*COEF(ISUB,5)/(TAU+TAUR2)+ 
     &    (ATAU0/ATAU5)*COEF(ISUB,6)*TAU/((TAU-TAUR2)**2+GAMR2**2)  
        ENDIF   
        COMFAC=COMFAC*ATAU0/(TAU*H1)    
      ENDIF 
      IF(MINT(43).EQ.4.AND.ISET(ISUB).NE.5) THEN    
        AYST0=YSTMAX-YSTMIN 
        AYST1=0.5*(YSTMAX-YSTMIN)**2    
        AYST2=AYST1 
        AYST3=2.*(ATAN(EXP(YSTMAX))-ATAN(EXP(YSTMIN)))  
        H2=(AYST0/AYST1)*COEF(ISUB,7)*(YST-YSTMIN)+(AYST0/AYST2)*   
     &  COEF(ISUB,8)*(YSTMAX-YST)+(AYST0/AYST3)*COEF(ISUB,9)/COSH(YST)  
        COMFAC=COMFAC*AYST0/H2  
      ENDIF 
    
C...2 -> 1 processes: reduction in angular part of phase space integral 
C...for case of decaying resonance. 
      ACTH0=CTNMAX-CTNMIN+CTPMAX-CTPMIN 
      IF((ISET(ISUB).EQ.1.OR.ISET(ISUB).EQ.3)) THEN
         IF( KFPR(ISUB,1).GT.0) THEN
            IF (MDCY(KFPR(ISUB,1),1).EQ.1) THEN   
               IF(KFPR(ISUB,1).EQ.25.OR.KFPR(ISUB,1).EQ.37) THEN   
                  COMFAC=COMFAC*0.5*ACTH0   
               ELSE    
                  COMFAC=COMFAC*0.125*(3.*ACTH0+CTNMAX**3-CTNMIN**3+    
     &                 CTPMAX**3-CTPMIN**3)  
               ENDIF
            ENDIF
         ENDIF   
    
C...2 -> 2 processes: angular part of phase space integral. 
      ELSEIF(ISET(ISUB).EQ.2.OR.ISET(ISUB).EQ.4) THEN   
        ACTH1=LOG((MAX(RM34,RSQM-CTNMIN)*MAX(RM34,RSQM-CTPMIN))/    
     &  (MAX(RM34,RSQM-CTNMAX)*MAX(RM34,RSQM-CTPMAX)))  
        ACTH2=LOG((MAX(RM34,RSQM+CTNMAX)*MAX(RM34,RSQM+CTPMAX))/    
     &  (MAX(RM34,RSQM+CTNMIN)*MAX(RM34,RSQM+CTPMIN)))  
        ACTH3=1./MAX(RM34,RSQM-CTNMAX)-1./MAX(RM34,RSQM-CTNMIN)+    
     &  1./MAX(RM34,RSQM-CTPMAX)-1./MAX(RM34,RSQM-CTPMIN)   
        ACTH4=1./MAX(RM34,RSQM+CTNMIN)-1./MAX(RM34,RSQM+CTNMAX)+    
     &  1./MAX(RM34,RSQM+CTPMIN)-1./MAX(RM34,RSQM+CTPMAX)   
        H3=COEF(ISUB,10)+   
     &  (ACTH0/ACTH1)*COEF(ISUB,11)/MAX(RM34,RSQM-CTH)+ 
     &  (ACTH0/ACTH2)*COEF(ISUB,12)/MAX(RM34,RSQM+CTH)+ 
     &  (ACTH0/ACTH3)*COEF(ISUB,13)/MAX(RM34,RSQM-CTH)**2+  
     &  (ACTH0/ACTH4)*COEF(ISUB,14)/MAX(RM34,RSQM+CTH)**2   
        COMFAC=COMFAC*ACTH0*0.5*BE34/H3 
      ENDIF 
    
C...2 -> 3, 4 processes: phace space integral in tau'.  
      IF(MINT(43).GE.2.AND.(ISET(ISUB).EQ.3.OR.ISET(ISUB).EQ.4)) THEN   
        ATAUP0=LOG(TAUPMX/TAUPMN)   
        ATAUP1=((1.-TAU/TAUPMX)**4-(1.-TAU/TAUPMN)**4)/(4.*TAU) 
        H4=COEF(ISUB,15)+   
     &  ATAUP0/ATAUP1*COEF(ISUB,16)/TAUP*(1.-TAU/TAUP)**3   
        IF(1.-TAU/TAUP.GT.1.E-4) THEN   
          FZW=(1.+TAU/TAUP)*LOG(TAUP/TAU)-2.*(1.-TAU/TAUP)  
        ELSE    
          FZW=1./6.*(1.-TAU/TAUP)**3*TAU/TAUP   
        ENDIF   
        COMFAC=COMFAC*ATAUP0*FZW/H4 
      ENDIF 
    
C...Phase space integral for low-pT and multiple interactions.  
      IF(ISET(ISUB).EQ.5) THEN  
        COMFAC=PARU(1)*PARU(5)*FACK*0.5*VINT(2)/SH2 
        ATAU0=LOG(2.*(1.+SQRT(1.-XT2))/XT2-1.)  
        ATAU1=2.*ATAN(1./XT2-1.)/SQRT(XT2)  
        H1=COEF(ISUB,1)+(ATAU0/ATAU1)*COEF(ISUB,2)/SQRT(TAU)    
        COMFAC=COMFAC*ATAU0/H1  
        AYST0=YSTMAX-YSTMIN 
        AYST1=0.5*(YSTMAX-YSTMIN)**2    
        AYST3=2.*(ATAN(EXP(YSTMAX))-ATAN(EXP(YSTMIN)))  
        H2=(AYST0/AYST1)*COEF(ISUB,7)*(YST-YSTMIN)+(AYST0/AYST1)*   
     &  COEF(ISUB,8)*(YSTMAX-YST)+(AYST0/AYST3)*COEF(ISUB,9)/COSH(YST)  
        COMFAC=COMFAC*AYST0/H2  
        IF(MSTP(82).LE.1) COMFAC=COMFAC*XT2**2*(1./VINT(149)-1.)    
C...For MSTP(82)>=2 an additional factor (xT2/(xT2+VINT(149))**2 is 
C...introduced to make cross-section finite for xT2 -> 0.   
        IF(MSTP(82).GE.2) COMFAC=COMFAC*XT2**2/(VINT(149)*  
     &  (1.+VINT(149))) 
      ENDIF 
    
C...A: 2 -> 1, tree diagrams.   
    
  145 IF(ISUB.LE.10) THEN   
      IF(ISUB.EQ.1) THEN    
C...f + fb -> gamma*/Z0.    
        MINT(61)=2  
        CALL PYWIDT_HIJING(23,SQRT(SH),WDTP,WDTE)  
        FACZ=COMFAC*AEM**2*4./3.    
        DO 150 I=MINA,MAXA  
        IF(I.EQ.0.OR.KFAC(1,I)*KFAC(2,-I).EQ.0) GOTO 150    
        EI=KCHG(IABS(I),1)/3.   
        AI=SIGN(1.,EI)  
        VI=AI-4.*EI*XW  
        FACF=1. 
        IF(IABS(I).LE.10) FACF=FACA/3.  
        NCHN=NCHN+1 
        ISIG(NCHN,1)=I  
        ISIG(NCHN,2)=-I 
        ISIG(NCHN,3)=1  
        SIGH(NCHN)=FACF*FACZ*(EI**2*VINT(111)+EI*VI/(8.*XW*(1.-XW))*    
     &  SH*(SH-SQMZ)/((SH-SQMZ)**2+GMMZ**2)*VINT(112)+(VI**2+AI**2)/    
     &  (16.*XW*(1.-XW))**2*SH2/((SH-SQMZ)**2+GMMZ**2)*VINT(114))   
  150   CONTINUE    
    
      ELSEIF(ISUB.EQ.2) THEN    
C...f + fb' -> W+/-.    
        CALL PYWIDT_HIJING(24,SQRT(SH),WDTP,WDTE)  
        FACW=COMFAC*(AEM/XW)**2*1./24*SH2/((SH-SQMW)**2+GMMW**2)    
        DO 170 I=MIN1,MAX1  
        IF(I.EQ.0.OR.KFAC(1,I).EQ.0) GOTO 170   
        IA=IABS(I)  
        DO 160 J=MIN2,MAX2  
        IF(J.EQ.0.OR.KFAC(2,J).EQ.0) GOTO 160   
        JA=IABS(J)  
        IF(I*J.GT.0.OR.MOD(IA+JA,2).EQ.0) GOTO 160  
        IF((IA.LE.10.AND.JA.GT.10).OR.(IA.GT.10.AND.JA.LE.10)) GOTO 160 
        KCHW=(KCHG(IA,1)*ISIGN(1,I)+KCHG(JA,1)*ISIGN(1,J))/3    
        FACF=1. 
        IF(IA.LE.10) FACF=VCKM((IA+1)/2,(JA+1)/2)*FACA/3.   
        NCHN=NCHN+1 
        ISIG(NCHN,1)=I  
        ISIG(NCHN,2)=J  
        ISIG(NCHN,3)=1  
        SIGH(NCHN)=FACF*FACW*(WDTE(0,1)+WDTE(0,(5-KCHW)/2)+WDTE(0,4))   
  160   CONTINUE    
  170   CONTINUE    
    
      ELSEIF(ISUB.EQ.3) THEN    
C...f + fb -> H0.   
        CALL PYWIDT_HIJING(25,SQRT(SH),WDTP,WDTE)  
        FACH=COMFAC*(AEM/XW)**2*1./48.*(SH/SQMW)**2*    
     &  SH2/((SH-SQMH)**2+GMMH**2)*(WDTE(0,1)+WDTE(0,2)+WDTE(0,4))  
        DO 180 I=MINA,MAXA  
        IF(I.EQ.0.OR.KFAC(1,I)*KFAC(2,-I).EQ.0) GOTO 180    
        RMQ=PMAS(IABS(I),1)**2/SH   
        NCHN=NCHN+1 
        ISIG(NCHN,1)=I  
        ISIG(NCHN,2)=-I 
        ISIG(NCHN,3)=1  
        SIGH(NCHN)=FACH*RMQ*SQRT(MAX(0.,1.-4.*RMQ)) 
  180   CONTINUE    
    
      ELSEIF(ISUB.EQ.4) THEN    
C...gamma + W+/- -> W+/-.   
    
      ELSEIF(ISUB.EQ.5) THEN    
C...Z0 + Z0 -> H0.  
        CALL PYWIDT_HIJING(25,SQRT(SH),WDTP,WDTE)  
        FACH=COMFAC*1./(128.*PARU(1)**2*16.*(1.-XW)**3)*(AEM/XW)**4*    
     &  (SH/SQMW)**2*SH2/((SH-SQMH)**2+GMMH**2)*    
     &  (WDTE(0,1)+WDTE(0,2)+WDTE(0,4)) 
        DO 200 I=MIN1,MAX1  
        IF(I.EQ.0.OR.KFAC(1,I).EQ.0) GOTO 200   
        DO 190 J=MIN2,MAX2  
        IF(J.EQ.0.OR.KFAC(2,J).EQ.0) GOTO 190   
        EI=KCHG(IABS(I),1)/3.   
        AI=SIGN(1.,EI)  
        VI=AI-4.*EI*XW  
        EJ=KCHG(IABS(J),1)/3.   
        AJ=SIGN(1.,EJ)  
        VJ=AJ-4.*EJ*XW  
        NCHN=NCHN+1 
        ISIG(NCHN,1)=I  
        ISIG(NCHN,2)=J  
        ISIG(NCHN,3)=1  
        SIGH(NCHN)=FACH*(VI**2+AI**2)*(VJ**2+AJ**2) 
  190   CONTINUE    
  200   CONTINUE    
    
      ELSEIF(ISUB.EQ.6) THEN    
C...Z0 + W+/- -> W+/-.  
    
      ELSEIF(ISUB.EQ.7) THEN    
C...W+ + W- -> Z0.  
    
      ELSEIF(ISUB.EQ.8) THEN    
C...W+ + W- -> H0.  
        CALL PYWIDT_HIJING(25,SQRT(SH),WDTP,WDTE)  
        FACH=COMFAC*1./(128*PARU(1)**2)*(AEM/XW)**4*(SH/SQMW)**2*   
     &  SH2/((SH-SQMH)**2+GMMH**2)*(WDTE(0,1)+WDTE(0,2)+WDTE(0,4))  
        DO 220 I=MIN1,MAX1  
        IF(I.EQ.0.OR.KFAC(1,I).EQ.0) GOTO 220   
        EI=SIGN(1.,FLOAT(I))*KCHG(IABS(I),1)    
        DO 210 J=MIN2,MAX2  
        IF(J.EQ.0.OR.KFAC(2,J).EQ.0) GOTO 210   
        EJ=SIGN(1.,FLOAT(J))*KCHG(IABS(J),1)    
        IF(EI*EJ.GT.0.) GOTO 210    
        NCHN=NCHN+1 
        ISIG(NCHN,1)=I  
        ISIG(NCHN,2)=J  
        ISIG(NCHN,3)=1  
        SIGH(NCHN)=FACH*VINT(180+I)*VINT(180+J) 
  210   CONTINUE    
  220   CONTINUE    
      ENDIF 
    
C...B: 2 -> 2, tree diagrams.   
    
      ELSEIF(ISUB.LE.20) THEN   
      IF(ISUB.EQ.11) THEN   
C...f + f' -> f + f'.   
        FACQQ1=COMFAC*AS**2*4./9.*(SH2+UH2)/TH2 
        FACQQB=COMFAC*AS**2*4./9.*((SH2+UH2)/TH2*FACA-  
     &  MSTP(34)*2./3.*UH2/(SH*TH)) 
        FACQQ2=COMFAC*AS**2*4./9.*((SH2+TH2)/UH2-   
     &  MSTP(34)*2./3.*SH2/(TH*UH)) 
        DO 240 I=MIN1,MAX1  
        IF(I.EQ.0.OR.KFAC(1,I).EQ.0) GOTO 240   
        DO 230 J=MIN2,MAX2  
        IF(J.EQ.0.OR.KFAC(2,J).EQ.0) GOTO 230   
        NCHN=NCHN+1 
        ISIG(NCHN,1)=I  
        ISIG(NCHN,2)=J  
        ISIG(NCHN,3)=1  
        SIGH(NCHN)=FACQQ1   
        IF(I.EQ.-J) SIGH(NCHN)=FACQQB   
        IF(I.EQ.J) THEN 
          SIGH(NCHN)=0.5*SIGH(NCHN) 
          NCHN=NCHN+1   
          ISIG(NCHN,1)=I    
          ISIG(NCHN,2)=J    
          ISIG(NCHN,3)=2    
          SIGH(NCHN)=0.5*FACQQ2 
        ENDIF   
  230   CONTINUE    
  240   CONTINUE    
    
      ELSEIF(ISUB.EQ.12) THEN   
C...f + fb -> f' + fb' (q + qb -> q' + qb' only).   
        CALL PYWIDT_HIJING(21,SQRT(SH),WDTP,WDTE)  
        FACQQB=COMFAC*AS**2*4./9.*(TH2+UH2)/SH2*(WDTE(0,1)+WDTE(0,2)+   
     &  WDTE(0,3)+WDTE(0,4))    
        DO 250 I=MINA,MAXA  
        IF(I.EQ.0.OR.KFAC(1,I)*KFAC(2,-I).EQ.0) GOTO 250    
        NCHN=NCHN+1 
        ISIG(NCHN,1)=I  
        ISIG(NCHN,2)=-I 
        ISIG(NCHN,3)=1  
        SIGH(NCHN)=FACQQB   
  250   CONTINUE    
    
      ELSEIF(ISUB.EQ.13) THEN   
C...f + fb -> g + g (q + qb -> g + g only). 
        FACGG1=COMFAC*AS**2*32./27.*(UH/TH-(2.+MSTP(34)*1./4.)*UH2/SH2) 
        FACGG2=COMFAC*AS**2*32./27.*(TH/UH-(2.+MSTP(34)*1./4.)*TH2/SH2) 
        DO 260 I=MINA,MAXA  
        IF(I.EQ.0.OR.KFAC(1,I)*KFAC(2,-I).EQ.0) GOTO 260    
        NCHN=NCHN+1 
        ISIG(NCHN,1)=I  
        ISIG(NCHN,2)=-I 
        ISIG(NCHN,3)=1  
        SIGH(NCHN)=0.5*FACGG1   
        NCHN=NCHN+1 
        ISIG(NCHN,1)=I  
        ISIG(NCHN,2)=-I 
        ISIG(NCHN,3)=2  
        SIGH(NCHN)=0.5*FACGG2   
  260   CONTINUE    
    
      ELSEIF(ISUB.EQ.14) THEN   
C...f + fb -> g + gamma (q + qb -> g + gamma only). 
        FACGG=COMFAC*AS*AEM*8./9.*(TH2+UH2)/(TH*UH) 
        DO 270 I=MINA,MAXA  
        IF(I.EQ.0.OR.KFAC(1,I)*KFAC(2,-I).EQ.0) GOTO 270    
        EI=KCHG(IABS(I),1)/3.   
        NCHN=NCHN+1 
        ISIG(NCHN,1)=I  
        ISIG(NCHN,2)=-I 
        ISIG(NCHN,3)=1  
        SIGH(NCHN)=FACGG*EI**2  
  270   CONTINUE    
    
      ELSEIF(ISUB.EQ.15) THEN   
C...f + fb -> g + Z0 (q + qb -> g + Z0 only).   
        FACZG=COMFAC*AS*AEM/(XW*(1.-XW))*1./18.*    
     &  (TH2+UH2+2.*SQM4*SH)/(TH*UH)    
        FACZG=FACZG*WIDS(23,2)  
        DO 280 I=MINA,MAXA  
        IF(I.EQ.0.OR.KFAC(1,I)*KFAC(2,-I).EQ.0) GOTO 280    
        EI=KCHG(IABS(I),1)/3.   
        AI=SIGN(1.,EI)  
        VI=AI-4.*EI*XW  
        NCHN=NCHN+1 
        ISIG(NCHN,1)=I  
        ISIG(NCHN,2)=-I 
        ISIG(NCHN,3)=1  
        SIGH(NCHN)=FACZG*(VI**2+AI**2)  
  280   CONTINUE    
    
      ELSEIF(ISUB.EQ.16) THEN   
C...f + fb' -> g + W+/- (q + qb' -> g + W+/- only). 
        FACWG=COMFAC*AS*AEM/XW*2./9.*(TH2+UH2+2.*SQM4*SH)/(TH*UH)   
        DO 300 I=MIN1,MAX1  
        IF(I.EQ.0.OR.KFAC(1,I).EQ.0) GOTO 300   
        IA=IABS(I)  
        DO 290 J=MIN2,MAX2  
        IF(J.EQ.0.OR.KFAC(2,J).EQ.0) GOTO 290   
        JA=IABS(J)  
        IF(I*J.GT.0.OR.MOD(IA+JA,2).EQ.0) GOTO 290  
        KCHW=(KCHG(IA,1)*ISIGN(1,I)+KCHG(JA,1)*ISIGN(1,J))/3    
        FCKM=1. 
        IF(MINT(43).EQ.4) FCKM=VCKM((IA+1)/2,(JA+1)/2)  
        NCHN=NCHN+1 
        ISIG(NCHN,1)=I  
        ISIG(NCHN,2)=J  
        ISIG(NCHN,3)=1  
        SIGH(NCHN)=FACWG*FCKM*WIDS(24,(5-KCHW)/2)   
  290   CONTINUE    
  300   CONTINUE    
    
      ELSEIF(ISUB.EQ.17) THEN   
C...f + fb -> g + H0 (q + qb -> g + H0 only).   
    
      ELSEIF(ISUB.EQ.18) THEN   
C...f + fb -> gamma + gamma.    
        FACGG=COMFAC*FACA*AEM**2*1./3.*(TH2+UH2)/(TH*UH)    
        DO 310 I=MINA,MAXA  
        IF(I.EQ.0.OR.KFAC(1,I)*KFAC(2,-I).EQ.0) GOTO 310    
        EI=KCHG(IABS(I),1)/3.   
        NCHN=NCHN+1 
        ISIG(NCHN,1)=I  
        ISIG(NCHN,2)=-I 
        ISIG(NCHN,3)=1  
        SIGH(NCHN)=FACGG*EI**4  
  310   CONTINUE    
    
      ELSEIF(ISUB.EQ.19) THEN   
C...f + fb -> gamma + Z0.   
        FACGZ=COMFAC*FACA*AEM**2/(XW*(1.-XW))*1./24.*   
     &  (TH2+UH2+2.*SQM4*SH)/(TH*UH)    
        FACGZ=FACGZ*WIDS(23,2)  
        DO 320 I=MINA,MAXA  
        IF(I.EQ.0.OR.KFAC(1,I)*KFAC(2,-I).EQ.0) GOTO 320    
        EI=KCHG(IABS(I),1)/3.   
        AI=SIGN(1.,EI)  
        VI=AI-4.*EI*XW  
        NCHN=NCHN+1 
        ISIG(NCHN,1)=I  
        ISIG(NCHN,2)=-I 
        ISIG(NCHN,3)=1  
        SIGH(NCHN)=FACGZ*EI**2*(VI**2+AI**2)    
  320   CONTINUE    
    
      ELSEIF(ISUB.EQ.20) THEN   
C...f + fb' -> gamma + W+/-.    
        FACGW=COMFAC*FACA*AEM**2/XW*1./6.*  
     &  ((2.*UH-TH)/(3.*(SH-SQM4)))**2*(TH2+UH2+2.*SQM4*SH)/(TH*UH) 
        DO 340 I=MIN1,MAX1  
        IF(I.EQ.0.OR.KFAC(1,I).EQ.0) GOTO 340   
        IA=IABS(I)  
        DO 330 J=MIN2,MAX2  
        IF(J.EQ.0.OR.KFAC(2,J).EQ.0) GOTO 330   
        JA=IABS(J)  
        IF(I*J.GT.0.OR.MOD(IA+JA,2).EQ.0) GOTO 330  
        KCHW=(KCHG(IA,1)*ISIGN(1,I)+KCHG(JA,1)*ISIGN(1,J))/3    
        FCKM=1. 
        IF(MINT(43).EQ.4) FCKM=VCKM((IA+1)/2,(JA+1)/2)  
        NCHN=NCHN+1 
        ISIG(NCHN,1)=I  
        ISIG(NCHN,2)=J  
        ISIG(NCHN,3)=1  
        SIGH(NCHN)=FACGW*FCKM*WIDS(24,(5-KCHW)/2)   
  330   CONTINUE    
  340   CONTINUE    
      ENDIF 
    
      ELSEIF(ISUB.LE.30) THEN   
      IF(ISUB.EQ.21) THEN   
C...f + fb -> gamma + H0.   
    
      ELSEIF(ISUB.EQ.22) THEN   
C...f + fb -> Z0 + Z0.  
        FACZZ=COMFAC*FACA*(AEM/(XW*(1.-XW)))**2*1./768.*    
     &  (UH/TH+TH/UH+2.*(SQM3+SQM4)*SH/(TH*UH)- 
     &  SQM3*SQM4*(1./TH2+1./UH2))  
        FACZZ=FACZZ*WIDS(23,1)  
        DO 350 I=MINA,MAXA  
        IF(I.EQ.0.OR.KFAC(1,I)*KFAC(2,-I).EQ.0) GOTO 350    
        EI=KCHG(IABS(I),1)/3.   
        AI=SIGN(1.,EI)  
        VI=AI-4.*EI*XW  
        NCHN=NCHN+1 
        ISIG(NCHN,1)=I  
        ISIG(NCHN,2)=-I 
        ISIG(NCHN,3)=1  
        SIGH(NCHN)=FACZZ*(VI**4+6.*VI**2*AI**2+AI**4)   
  350   CONTINUE    
    
      ELSEIF(ISUB.EQ.23) THEN   
C...f + fb' -> Z0 + W+/-.   
        FACZW=COMFAC*FACA*(AEM/XW)**2*1./6. 
        FACZW=FACZW*WIDS(23,2)  
        THUH=MAX(TH*UH-SQM3*SQM4,SH*CKIN(3)**2) 
        DO 370 I=MIN1,MAX1  
        IF(I.EQ.0.OR.KFAC(1,I).EQ.0) GOTO 370   
        IA=IABS(I)  
        DO 360 J=MIN2,MAX2  
        IF(J.EQ.0.OR.KFAC(2,J).EQ.0) GOTO 360   
        JA=IABS(J)  
        IF(I*J.GT.0.OR.MOD(IA+JA,2).EQ.0) GOTO 360  
        KCHW=(KCHG(IA,1)*ISIGN(1,I)+KCHG(JA,1)*ISIGN(1,J))/3    
        EI=KCHG(IA,1)/3.    
        AI=SIGN(1.,EI)  
        VI=AI-4.*EI*XW  
        EJ=KCHG(JA,1)/3.    
        AJ=SIGN(1.,EJ)  
        VJ=AJ-4.*EJ*XW  
        IF(VI+AI.GT.0) THEN 
          VISAV=VI  
          AISAV=AI  
          VI=VJ 
          AI=AJ 
          VJ=VISAV  
          AJ=AISAV  
        ENDIF   
        FCKM=1. 
        IF(MINT(43).EQ.4) FCKM=VCKM((IA+1)/2,(JA+1)/2)  
        NCHN=NCHN+1 
        ISIG(NCHN,1)=I  
        ISIG(NCHN,2)=J  
        ISIG(NCHN,3)=1  
        SIGH(NCHN)=FACZW*FCKM*(1./(SH-SQMW)**2* 
     &  ((9.-8.*XW)/4.*THUH+(8.*XW-6.)/4.*SH*(SQM3+SQM4))+  
     &  (THUH-SH*(SQM3+SQM4))/(2.*(SH-SQMW))*((VJ+AJ)/TH-(VI+AI)/UH)+   
     &  THUH/(16.*(1.-XW))*((VJ+AJ)**2/TH2+(VI+AI)**2/UH2)+ 
     &  SH*(SQM3+SQM4)/(8.*(1.-XW))*(VI+AI)*(VJ+AJ)/(TH*UH))*   
     &  WIDS(24,(5-KCHW)/2) 
  360   CONTINUE    
  370   CONTINUE    
    
      ELSEIF(ISUB.EQ.24) THEN   
C...f + fb -> Z0 + H0.  
        THUH=MAX(TH*UH-SQM3*SQM4,SH*CKIN(3)**2) 
        FACHZ=COMFAC*FACA*(AEM/(XW*(1.-XW)))**2*1./96.* 
     &  (THUH+2.*SH*SQMZ)/(SH-SQMZ)**2  
        FACHZ=FACHZ*WIDS(23,2)*WIDS(25,2)   
        DO 380 I=MINA,MAXA  
        IF(I.EQ.0.OR.KFAC(1,I)*KFAC(2,-I).EQ.0) GOTO 380    
        EI=KCHG(IABS(I),1)/3.   
        AI=SIGN(1.,EI)  
        VI=AI-4.*EI*XW  
        NCHN=NCHN+1 
        ISIG(NCHN,1)=I  
        ISIG(NCHN,2)=-I 
        ISIG(NCHN,3)=1  
        SIGH(NCHN)=FACHZ*(VI**2+AI**2)  
  380   CONTINUE    
    
      ELSEIF(ISUB.EQ.25) THEN   
C...f + fb -> W+ + W-.  
        FACWW=COMFAC*FACA*(AEM/XW)**2*1./12.    
        FACWW=FACWW*WIDS(24,1)  
        THUH=MAX(TH*UH-SQM3*SQM4,SH*CKIN(3)**2) 
        DO 390 I=MINA,MAXA  
        IF(I.EQ.0.OR.KFAC(1,I)*KFAC(2,-I).EQ.0) GOTO 390    
        EI=KCHG(IABS(I),1)/3.   
        AI=SIGN(1.,EI)  
        VI=AI-4.*EI*XW  
        DSIGWW=THUH/SH2*(3.-(SH-3.*(SQM3+SQM4))/(SH-SQMZ)*  
     &  (VI+AI)/(2.*AI*(1.-XW))+(SH/(SH-SQMZ))**2*  
     &  (1.-2.*(SQM3+SQM4)/SH+12.*SQM3*SQM4/SH2)*(VI**2+AI**2)/ 
     &  (8.*(1.-XW)**2))-2.*SQMZ/(SH-SQMZ)*(VI+AI)/AI+  
     &  SQMZ*SH/(SH-SQMZ)**2*(1.-2.*(SQM3+SQM4)/SH)*(VI**2+AI**2)/  
     &  (2.*(1.-XW))    
        IF(KCHG(IABS(I),1).LT.0) THEN   
          DSIGWW=DSIGWW+2.*(1.+SQMZ/(SH-SQMZ)*(VI+AI)/(2.*AI))* 
     &    (THUH/(SH*TH)-(SQM3+SQM4)/TH)+THUH/TH2    
        ELSE    
          DSIGWW=DSIGWW+2.*(1.+SQMZ/(SH-SQMZ)*(VI+AI)/(2.*AI))* 
     &    (THUH/(SH*UH)-(SQM3+SQM4)/UH)+THUH/UH2    
        ENDIF   
        NCHN=NCHN+1 
        ISIG(NCHN,1)=I  
        ISIG(NCHN,2)=-I 
        ISIG(NCHN,3)=1  
        SIGH(NCHN)=FACWW*DSIGWW 
  390   CONTINUE    
    
      ELSEIF(ISUB.EQ.26) THEN   
C...f + fb' -> W+/- + H0.   
        THUH=MAX(TH*UH-SQM3*SQM4,SH*CKIN(3)**2) 
        FACHW=COMFAC*FACA*(AEM/XW)**2*1./24.*(THUH+2.*SH*SQMW)/ 
     &  (SH-SQMW)**2    
        FACHW=FACHW*WIDS(25,2)  
        DO 410 I=MIN1,MAX1  
        IF(I.EQ.0.OR.KFAC(1,I).EQ.0) GOTO 410   
        IA=IABS(I)  
        DO 400 J=MIN2,MAX2  
        IF(J.EQ.0.OR.KFAC(1,J).EQ.0) GOTO 400   
        JA=IABS(J)  
        IF(I*J.GT.0.OR.MOD(IA+JA,2).EQ.0) GOTO 400  
        KCHW=(KCHG(IA,1)*ISIGN(1,I)+KCHG(JA,1)*ISIGN(1,J))/3    
        FCKM=1. 
        IF(MINT(43).EQ.4) FCKM=VCKM((IA+1)/2,(JA+1)/2)  
        NCHN=NCHN+1 
        ISIG(NCHN,1)=I  
        ISIG(NCHN,2)=J  
        ISIG(NCHN,3)=1  
        SIGH(NCHN)=FACHW*FCKM*WIDS(24,(5-KCHW)/2)   
  400   CONTINUE    
  410   CONTINUE    
    
      ELSEIF(ISUB.EQ.27) THEN   
C...f + fb -> H0 + H0.  
    
      ELSEIF(ISUB.EQ.28) THEN   
C...f + g -> f + g (q + g -> q + g only).   
        FACQG1=COMFAC*AS**2*4./9.*((2.+MSTP(34)*1./4.)*UH2/TH2-UH/SH)*  
     &  FACA    
        FACQG2=COMFAC*AS**2*4./9.*((2.+MSTP(34)*1./4.)*SH2/TH2-SH/UH)   
        DO 430 I=MINA,MAXA  
        IF(I.EQ.0) GOTO 430 
        DO 420 ISDE=1,2 
        IF(ISDE.EQ.1.AND.KFAC(1,I)*KFAC(2,21).EQ.0) GOTO 420    
        IF(ISDE.EQ.2.AND.KFAC(1,21)*KFAC(2,I).EQ.0) GOTO 420    
        NCHN=NCHN+1 
        ISIG(NCHN,ISDE)=I   
        ISIG(NCHN,3-ISDE)=21    
        ISIG(NCHN,3)=1  
        SIGH(NCHN)=FACQG1   
        NCHN=NCHN+1 
        ISIG(NCHN,ISDE)=I   
        ISIG(NCHN,3-ISDE)=21    
        ISIG(NCHN,3)=2  
        SIGH(NCHN)=FACQG2   
  420   CONTINUE    
  430   CONTINUE    
    
      ELSEIF(ISUB.EQ.29) THEN   
C...f + g -> f + gamma (q + g -> q + gamma only).   
        FGQ=COMFAC*FACA*AS*AEM*1./3.*(SH2+UH2)/(-SH*UH) 
        DO 450 I=MINA,MAXA  
        IF(I.EQ.0) GOTO 450 
        EI=KCHG(IABS(I),1)/3.   
        FACGQ=FGQ*EI**2 
        DO 440 ISDE=1,2 
        IF(ISDE.EQ.1.AND.KFAC(1,I)*KFAC(2,21).EQ.0) GOTO 440    
        IF(ISDE.EQ.2.AND.KFAC(1,21)*KFAC(2,I).EQ.0) GOTO 440    
        NCHN=NCHN+1 
        ISIG(NCHN,ISDE)=I   
        ISIG(NCHN,3-ISDE)=21    
        ISIG(NCHN,3)=1  
        SIGH(NCHN)=FACGQ    
  440   CONTINUE    
  450   CONTINUE    
    
      ELSEIF(ISUB.EQ.30) THEN   
C...f + g -> f + Z0 (q + g -> q + Z0 only). 
        FZQ=COMFAC*FACA*AS*AEM/(XW*(1.-XW))*1./48.* 
     &  (SH2+UH2+2.*SQM4*TH)/(-SH*UH)   
        FZQ=FZQ*WIDS(23,2)  
        DO 470 I=MINA,MAXA  
        IF(I.EQ.0) GOTO 470 
        EI=KCHG(IABS(I),1)/3.   
        AI=SIGN(1.,EI)  
        VI=AI-4.*EI*XW  
        FACZQ=FZQ*(VI**2+AI**2) 
        DO 460 ISDE=1,2 
        IF(ISDE.EQ.1.AND.KFAC(1,I)*KFAC(2,21).EQ.0) GOTO 460    
        IF(ISDE.EQ.2.AND.KFAC(1,21)*KFAC(2,I).EQ.0) GOTO 460    
        NCHN=NCHN+1 
        ISIG(NCHN,ISDE)=I   
        ISIG(NCHN,3-ISDE)=21    
        ISIG(NCHN,3)=1  
        SIGH(NCHN)=FACZQ    
  460   CONTINUE    
  470   CONTINUE    
      ENDIF 
    
      ELSEIF(ISUB.LE.40) THEN   
      IF(ISUB.EQ.31) THEN   
C...f + g -> f' + W+/- (q + g -> q' + W+/- only).   
        FACWQ=COMFAC*FACA*AS*AEM/XW*1./12.* 
     &  (SH2+UH2+2.*SQM4*TH)/(-SH*UH)   
        DO 490 I=MINA,MAXA  
        IF(I.EQ.0) GOTO 490 
        IA=IABS(I)  
        KCHW=ISIGN(1,KCHG(IA,1)*ISIGN(1,I)) 
        DO 480 ISDE=1,2 
        IF(ISDE.EQ.1.AND.KFAC(1,I)*KFAC(2,21).EQ.0) GOTO 480    
        IF(ISDE.EQ.2.AND.KFAC(1,21)*KFAC(2,I).EQ.0) GOTO 480    
        NCHN=NCHN+1 
        ISIG(NCHN,ISDE)=I   
        ISIG(NCHN,3-ISDE)=21    
        ISIG(NCHN,3)=1  
        SIGH(NCHN)=FACWQ*VINT(180+I)*WIDS(24,(5-KCHW)/2)    
  480   CONTINUE    
  490   CONTINUE    
    
      ELSEIF(ISUB.EQ.32) THEN   
C...f + g -> f + H0 (q + g -> q + H0 only). 
    
      ELSEIF(ISUB.EQ.33) THEN   
C...f + gamma -> f + g (q + gamma -> q + g only).   
    
      ELSEIF(ISUB.EQ.34) THEN   
C...f + gamma -> f + gamma. 
    
      ELSEIF(ISUB.EQ.35) THEN   
C...f + gamma -> f + Z0.    
    
      ELSEIF(ISUB.EQ.36) THEN   
C...f + gamma -> f' + W+/-. 
    
      ELSEIF(ISUB.EQ.37) THEN   
C...f + gamma -> f + H0.    
    
      ELSEIF(ISUB.EQ.38) THEN   
C...f + Z0 -> f + g (q + Z0 -> q + g only). 
    
      ELSEIF(ISUB.EQ.39) THEN   
C...f + Z0 -> f + gamma.    
    
      ELSEIF(ISUB.EQ.40) THEN   
C...f + Z0 -> f + Z0.   
      ENDIF 
    
      ELSEIF(ISUB.LE.50) THEN   
      IF(ISUB.EQ.41) THEN   
C...f + Z0 -> f' + W+/-.    
    
      ELSEIF(ISUB.EQ.42) THEN   
C...f + Z0 -> f + H0.   
    
      ELSEIF(ISUB.EQ.43) THEN   
C...f + W+/- -> f' + g (q + W+/- -> q' + g only).   
    
      ELSEIF(ISUB.EQ.44) THEN   
C...f + W+/- -> f' + gamma. 
    
      ELSEIF(ISUB.EQ.45) THEN   
C...f + W+/- -> f' + Z0.    
    
      ELSEIF(ISUB.EQ.46) THEN   
C...f + W+/- -> f' + W+/-.  
    
      ELSEIF(ISUB.EQ.47) THEN   
C...f + W+/- -> f' + H0.    
    
      ELSEIF(ISUB.EQ.48) THEN   
C...f + H0 -> f + g (q + H0 -> q + g only). 
    
      ELSEIF(ISUB.EQ.49) THEN   
C...f + H0 -> f + gamma.    
    
      ELSEIF(ISUB.EQ.50) THEN   
C...f + H0 -> f + Z0.   
      ENDIF 
    
      ELSEIF(ISUB.LE.60) THEN   
      IF(ISUB.EQ.51) THEN   
C...f + H0 -> f' + W+/-.    
    
      ELSEIF(ISUB.EQ.52) THEN   
C...f + H0 -> f + H0.   
    
      ELSEIF(ISUB.EQ.53) THEN   
C...g + g -> f + fb (g + g -> q + qb only). 
        CALL PYWIDT_HIJING(21,SQRT(SH),WDTP,WDTE)  
        FACQQ1=COMFAC*AS**2*1./6.*(UH/TH-(2.+MSTP(34)*1./4.)*UH2/SH2)*  
     &  (WDTE(0,1)+WDTE(0,2)+WDTE(0,3)+WDTE(0,4))*FACA  
        FACQQ2=COMFAC*AS**2*1./6.*(TH/UH-(2.+MSTP(34)*1./4.)*TH2/SH2)*  
     &  (WDTE(0,1)+WDTE(0,2)+WDTE(0,3)+WDTE(0,4))*FACA  
        IF(KFAC(1,21)*KFAC(2,21).EQ.0) GOTO 500 
        NCHN=NCHN+1 
        ISIG(NCHN,1)=21 
        ISIG(NCHN,2)=21 
        ISIG(NCHN,3)=1  
        SIGH(NCHN)=FACQQ1   
        NCHN=NCHN+1 
        ISIG(NCHN,1)=21 
        ISIG(NCHN,2)=21 
        ISIG(NCHN,3)=2  
        SIGH(NCHN)=FACQQ2   
  500   CONTINUE    
    
      ELSEIF(ISUB.EQ.54) THEN   
C...g + gamma -> f + fb (g + gamma -> q + qb only). 
    
      ELSEIF(ISUB.EQ.55) THEN   
C...g + gamma -> f + fb (g + gamma -> q + qb only). 
    
      ELSEIF(ISUB.EQ.56) THEN   
C...g + gamma -> f + fb (g + gamma -> q + qb only). 
    
      ELSEIF(ISUB.EQ.57) THEN   
C...g + gamma -> f + fb (g + gamma -> q + qb only). 
    
      ELSEIF(ISUB.EQ.58) THEN   
C...gamma + gamma -> f + fb.    
    
      ELSEIF(ISUB.EQ.59) THEN   
C...gamma + Z0 -> f + fb.   
    
      ELSEIF(ISUB.EQ.60) THEN   
C...gamma + W+/- -> f + fb'.    
      ENDIF 
    
      ELSEIF(ISUB.LE.70) THEN   
      IF(ISUB.EQ.61) THEN   
C...gamma + H0 -> f + fb.   
    
      ELSEIF(ISUB.EQ.62) THEN   
C...Z0 + Z0 -> f + fb.  
    
      ELSEIF(ISUB.EQ.63) THEN   
C...Z0 + W+/- -> f + fb'.   
    
      ELSEIF(ISUB.EQ.64) THEN   
C...Z0 + H0 -> f + fb.  
    
      ELSEIF(ISUB.EQ.65) THEN   
C...W+ + W- -> f + fb.  
    
      ELSEIF(ISUB.EQ.66) THEN   
C...W+/- + H0 -> f + fb'.   
    
      ELSEIF(ISUB.EQ.67) THEN   
C...H0 + H0 -> f + fb.  
    
      ELSEIF(ISUB.EQ.68) THEN   
C...g + g -> g + g. 
        FACGG1=COMFAC*AS**2*9./4.*(SH2/TH2+2.*SH/TH+3.+2.*TH/SH+    
     &  TH2/SH2)*FACA   
        FACGG2=COMFAC*AS**2*9./4.*(UH2/SH2+2.*UH/SH+3.+2.*SH/UH+    
     &  SH2/UH2)*FACA   
        FACGG3=COMFAC*AS**2*9./4.*(TH2/UH2+2.*TH/UH+3+2.*UH/TH+UH2/TH2) 
        IF(KFAC(1,21)*KFAC(2,21).EQ.0) GOTO 510 
        NCHN=NCHN+1 
        ISIG(NCHN,1)=21 
        ISIG(NCHN,2)=21 
        ISIG(NCHN,3)=1  
        SIGH(NCHN)=0.5*FACGG1   
        NCHN=NCHN+1 
        ISIG(NCHN,1)=21 
        ISIG(NCHN,2)=21 
        ISIG(NCHN,3)=2  
        SIGH(NCHN)=0.5*FACGG2   
        NCHN=NCHN+1 
        ISIG(NCHN,1)=21 
        ISIG(NCHN,2)=21 
        ISIG(NCHN,3)=3  
        SIGH(NCHN)=0.5*FACGG3   
  510   CONTINUE    
    
      ELSEIF(ISUB.EQ.69) THEN   
C...gamma + gamma -> W+ + W-.   
    
      ELSEIF(ISUB.EQ.70) THEN   
C...gamma + W+/- -> gamma + W+/-.   
      ENDIF 
    
      ELSEIF(ISUB.LE.80) THEN   
      IF(ISUB.EQ.71) THEN   
C...Z0 + Z0 -> Z0 + Z0. 
        BE2=1.-4.*SQMZ/SH   
        TH=-0.5*SH*BE2*(1.-CTH) 
        UH=-0.5*SH*BE2*(1.+CTH) 
        SHANG=1./(1.-XW)*SQMW/SQMZ*(1.+BE2)**2  
        ASHRE=(SH-SQMH)/((SH-SQMH)**2+GMMH**2)*SHANG    
        ASHIM=-GMMH/((SH-SQMH)**2+GMMH**2)*SHANG    
        THANG=1./(1.-XW)*SQMW/SQMZ*(BE2-CTH)**2 
        ATHRE=(TH-SQMH)/((TH-SQMH)**2+GMMH**2)*THANG    
        ATHIM=-GMMH/((TH-SQMH)**2+GMMH**2)*THANG    
        UHANG=1./(1.-XW)*SQMW/SQMZ*(BE2+CTH)**2 
        AUHRE=(UH-SQMH)/((UH-SQMH)**2+GMMH**2)*UHANG    
        AUHIM=-GMMH/((UH-SQMH)**2+GMMH**2)*UHANG    
        FACH=0.5*COMFAC*1./(4096.*PARU(1)**2*16.*(1.-XW)**2)*   
     &  (AEM/XW)**4*(SH/SQMW)**2*((ASHRE+ATHRE+AUHRE)**2+   
     &  (ASHIM+ATHIM+AUHIM)**2)*SQMZ/SQMW   
        DO 530 I=MIN1,MAX1  
        IF(I.EQ.0.OR.KFAC(1,I).EQ.0) GOTO 530   
        EI=KCHG(IABS(I),1)/3.   
        AI=SIGN(1.,EI)  
        VI=AI-4.*EI*XW  
        AVI=AI**2+VI**2 
        DO 520 J=MIN2,MAX2  
        IF(J.EQ.0.OR.KFAC(2,J).EQ.0) GOTO 520   
        EJ=KCHG(IABS(J),1)/3.   
        AJ=SIGN(1.,EJ)  
        VJ=AJ-4.*EJ*XW  
        AVJ=AJ**2+VJ**2 
        NCHN=NCHN+1 
        ISIG(NCHN,1)=I  
        ISIG(NCHN,2)=J  
        ISIG(NCHN,3)=1  
        SIGH(NCHN)=FACH*AVI*AVJ 
  520   CONTINUE    
  530   CONTINUE    
    
      ELSEIF(ISUB.EQ.72) THEN   
C...Z0 + Z0 -> W+ + W-. 
        BE2=SQRT((1.-4.*SQMW/SH)*(1.-4.*SQMZ/SH))   
        CTH2=CTH**2 
        TH=-0.5*SH*(1.-2.*(SQMW+SQMZ)/SH-BE2*CTH)   
        UH=-0.5*SH*(1.-2.*(SQMW+SQMZ)/SH+BE2*CTH)   
        SHANG=4.*SQRT(SQMW/(SQMZ*(1.-XW)))*(1.-2.*SQMW/SH)* 
     &  (1.-2.*SQMZ/SH) 
        ASHRE=(SH-SQMH)/((SH-SQMH)**2+GMMH**2)*SHANG    
        ASHIM=-GMMH/((SH-SQMH)**2+GMMH**2)*SHANG    
        ATWRE=(1.-XW)/SQMZ*SH/(TH-SQMW)*((CTH-BE2)**2*(3./2.+BE2/2.*CTH-    
     &  (SQMW+SQMZ)/SH+(SQMW-SQMZ)**2/(SH*SQMW))+4.*((SQMW+SQMZ)/SH*    
     &  (1.-3.*CTH2)+8.*SQMW*SQMZ/SH2*(2.*CTH2-1.)+ 
     &  4.*(SQMW**2+SQMZ**2)/SH2*CTH2+2.*(SQMW+SQMZ)/SH*BE2*CTH))   
        ATWIM=0.    
        AUWRE=(1.-XW)/SQMZ*SH/(UH-SQMW)*((CTH+BE2)**2*(3./2.-BE2/2.*CTH-    
     &  (SQMW+SQMZ)/SH+(SQMW-SQMZ)**2/(SH*SQMW))+4.*((SQMW+SQMZ)/SH*    
     &  (1.-3.*CTH2)+8.*SQMW*SQMZ/SH2*(2.*CTH2-1.)+ 
     &  4.*(SQMW**2+SQMZ**2)/SH2*CTH2-2.*(SQMW+SQMZ)/SH*BE2*CTH))   
        AUWIM=0.    
        A4RE=2.*(1.-XW)/SQMZ*(3.-CTH2-4.*(SQMW+SQMZ)/SH)    
        A4IM=0. 
        FACH=COMFAC*1./(4096.*PARU(1)**2*16.*(1.-XW)**2)*(AEM/XW)**4*   
     &  (SH/SQMW)**2*((ASHRE+ATWRE+AUWRE+A4RE)**2+  
     &  (ASHIM+ATWIM+AUWIM+A4IM)**2)*SQMZ/SQMW  
        DO 550 I=MIN1,MAX1  
        IF(I.EQ.0.OR.KFAC(1,I).EQ.0) GOTO 550   
        EI=KCHG(IABS(I),1)/3.   
        AI=SIGN(1.,EI)  
        VI=AI-4.*EI*XW  
        AVI=AI**2+VI**2 
        DO 540 J=MIN2,MAX2  
        IF(J.EQ.0.OR.KFAC(2,J).EQ.0) GOTO 540   
        EJ=KCHG(IABS(J),1)/3.   
        AJ=SIGN(1.,EJ)  
        VJ=AJ-4.*EJ*XW  
        AVJ=AJ**2+VJ**2 
        NCHN=NCHN+1 
        ISIG(NCHN,1)=I  
        ISIG(NCHN,2)=J  
        ISIG(NCHN,3)=1  
        SIGH(NCHN)=FACH*AVI*AVJ 
  540   CONTINUE    
  550   CONTINUE    
    
      ELSEIF(ISUB.EQ.73) THEN   
C...Z0 + W+/- -> Z0 + W+/-. 
        BE2=1.-2.*(SQMZ+SQMW)/SH+((SQMZ-SQMW)/SH)**2    
        EP1=1.+(SQMZ-SQMW)/SH   
        EP2=1.-(SQMZ-SQMW)/SH   
        TH=-0.5*SH*BE2*(1.-CTH) 
        UH=(SQMZ-SQMW)**2/SH-0.5*SH*BE2*(1.+CTH)    
        THANG=SQRT(SQMW/(SQMZ*(1.-XW)))*(BE2-EP1*CTH)*(BE2-EP2*CTH) 
        ATHRE=(TH-SQMH)/((TH-SQMH)**2+GMMH**2)*THANG    
        ATHIM=-GMMH/((TH-SQMH)**2+GMMH**2)*THANG    
        ASWRE=(1.-XW)/SQMZ*SH/(SH-SQMW)*(-BE2*(EP1+EP2)**4*CTH+ 
     &  1./4.*(BE2+EP1*EP2)**2*((EP1-EP2)**2-4.*BE2*CTH)+   
     &  2.*BE2*(BE2+EP1*EP2)*(EP1+EP2)**2*CTH-  
     &  1./16.*SH/SQMW*(EP1**2-EP2**2)**2*(BE2+EP1*EP2)**2) 
        ASWIM=0.    
        AUWRE=(1.-XW)/SQMZ*SH/(UH-SQMW)*(-BE2*(EP2+EP1*CTH)*    
     &  (EP1+EP2*CTH)*(BE2+EP1*EP2)+BE2*(EP2+EP1*CTH)*  
     &  (BE2+EP1*EP2*CTH)*(2.*EP2-EP2*CTH+EP1)-BE2*(EP2+EP1*CTH)**2*    
     &  (BE2-EP2**2*CTH)-1./8.*(BE2+EP1*EP2*CTH)**2*((EP1+EP2)**2+  
     &  2.*BE2*(1.-CTH))+1./32.*SH/SQMW*(BE2+EP1*EP2*CTH)**2*   
     &  (EP1**2-EP2**2)**2-BE2*(EP1+EP2*CTH)*(EP2+EP1*CTH)* 
     &  (BE2+EP1*EP2)+BE2*(EP1+EP2*CTH)*(BE2+EP1*EP2*CTH)*  
     &  (2.*EP1-EP1*CTH+EP2)-BE2*(EP1+EP2*CTH)**2*(BE2-EP1**2*CTH)- 
     &  1./8.*(BE2+EP1*EP2*CTH)**2*((EP1+EP2)**2+2.*BE2*(1.-CTH))+  
     &  1./32.*SH/SQMW*(BE2+EP1*EP2*CTH)**2*(EP1**2-EP2**2)**2) 
        AUWIM=0.    
        A4RE=(1.-XW)/SQMZ*(EP1**2*EP2**2*(CTH**2-1.)-   
     &  2.*BE2*(EP1**2+EP2**2+EP1*EP2)*CTH-2.*BE2*EP1*EP2)  
        A4IM=0. 
        FACH=COMFAC*1./(4096.*PARU(1)**2*4.*(1.-XW))*(AEM/XW)**4*   
     &  (SH/SQMW)**2*((ATHRE+ASWRE+AUWRE+A4RE)**2+  
     &  (ATHIM+ASWIM+AUWIM+A4IM)**2)*SQRT(SQMZ/SQMW)    
        DO 570 I=MIN1,MAX1  
        IF(I.EQ.0.OR.KFAC(1,I).EQ.0) GOTO 570   
        EI=KCHG(IABS(I),1)/3.   
        AI=SIGN(1.,EI)  
        VI=AI-4.*EI*XW  
        AVI=AI**2+VI**2 
        DO 560 J=MIN2,MAX2  
        IF(J.EQ.0.OR.KFAC(2,J).EQ.0) GOTO 560   
        EJ=KCHG(IABS(J),1)/3.   
        AJ=SIGN(1.,EJ)  
        VJ=AI-4.*EJ*XW  
        AVJ=AJ**2+VJ**2 
        NCHN=NCHN+1 
        ISIG(NCHN,1)=I  
        ISIG(NCHN,2)=J  
        ISIG(NCHN,3)=1  
        SIGH(NCHN)=FACH*(AVI*VINT(180+J)+VINT(180+I)*AVJ)   
  560   CONTINUE    
  570   CONTINUE    
    
      ELSEIF(ISUB.EQ.75) THEN   
C...W+ + W- -> gamma + gamma.   
    
      ELSEIF(ISUB.EQ.76) THEN   
C...W+ + W- -> Z0 + Z0. 
        BE2=SQRT((1.-4.*SQMW/SH)*(1.-4.*SQMZ/SH))   
        CTH2=CTH**2 
        TH=-0.5*SH*(1.-2.*(SQMW+SQMZ)/SH-BE2*CTH)   
        UH=-0.5*SH*(1.-2.*(SQMW+SQMZ)/SH+BE2*CTH)   
        SHANG=4.*SQRT(SQMW/(SQMZ*(1.-XW)))*(1.-2.*SQMW/SH)* 
     &  (1.-2.*SQMZ/SH) 
        ASHRE=(SH-SQMH)/((SH-SQMH)**2+GMMH**2)*SHANG    
        ASHIM=-GMMH/((SH-SQMH)**2+GMMH**2)*SHANG    
        ATWRE=(1.-XW)/SQMZ*SH/(TH-SQMW)*((CTH-BE2)**2*(3./2.+BE2/2.*CTH-    
     &  (SQMW+SQMZ)/SH+(SQMW-SQMZ)**2/(SH*SQMW))+4.*((SQMW+SQMZ)/SH*    
     &  (1.-3.*CTH2)+8.*SQMW*SQMZ/SH2*(2.*CTH2-1.)+ 
     &  4.*(SQMW**2+SQMZ**2)/SH2*CTH2+2.*(SQMW+SQMZ)/SH*BE2*CTH))   
        ATWIM=0.    
        AUWRE=(1.-XW)/SQMZ*SH/(UH-SQMW)*((CTH+BE2)**2*(3./2.-BE2/2.*CTH-    
     &  (SQMW+SQMZ)/SH+(SQMW-SQMZ)**2/(SH*SQMW))+4.*((SQMW+SQMZ)/SH*    
     &  (1.-3.*CTH2)+8.*SQMW*SQMZ/SH2*(2.*CTH2-1.)+ 
     &  4.*(SQMW**2+SQMZ**2)/SH2*CTH2-2.*(SQMW+SQMZ)/SH*BE2*CTH))   
        AUWIM=0.    
        A4RE=2.*(1.-XW)/SQMZ*(3.-CTH2-4.*(SQMW+SQMZ)/SH)    
        A4IM=0. 
        FACH=0.5*COMFAC*1./(4096.*PARU(1)**2)*(AEM/XW)**4*(SH/SQMW)**2* 
     &  ((ASHRE+ATWRE+AUWRE+A4RE)**2+(ASHIM+ATWIM+AUWIM+A4IM)**2)   
        DO 590 I=MIN1,MAX1  
        IF(I.EQ.0.OR.KFAC(1,I).EQ.0) GOTO 590   
        EI=SIGN(1.,FLOAT(I))*KCHG(IABS(I),1)    
        DO 580 J=MIN2,MAX2  
        IF(J.EQ.0.OR.KFAC(2,J).EQ.0) GOTO 580   
        EJ=SIGN(1.,FLOAT(J))*KCHG(IABS(J),1)    
        IF(EI*EJ.GT.0.) GOTO 580    
        NCHN=NCHN+1 
        ISIG(NCHN,1)=I  
        ISIG(NCHN,2)=J  
        ISIG(NCHN,3)=1  
        SIGH(NCHN)=FACH*VINT(180+I)*VINT(180+J) 
  580   CONTINUE    
  590   CONTINUE    
    
      ELSEIF(ISUB.EQ.77) THEN   
C...W+/- + W+/- -> W+/- + W+/-. 
        BE2=1.-4.*SQMW/SH   
        BE4=BE2**2  
        CTH2=CTH**2 
        CTH3=CTH**3 
        TH=-0.5*SH*BE2*(1.-CTH) 
        UH=-0.5*SH*BE2*(1.+CTH) 
        SHANG=(1.+BE2)**2   
        ASHRE=(SH-SQMH)/((SH-SQMH)**2+GMMH**2)*SHANG    
        ASHIM=-GMMH/((SH-SQMH)**2+GMMH**2)*SHANG    
        THANG=(BE2-CTH)**2  
        ATHRE=(TH-SQMH)/((TH-SQMH)**2+GMMH**2)*THANG    
        ATHIM=-GMMH/((TH-SQMH)**2+GMMH**2)*THANG    
        SGZANG=1./SQMW*BE2*(3.-BE2)**2*CTH  
        ASGRE=XW*SGZANG 
        ASGIM=0.    
        ASZRE=(1.-XW)*SH/(SH-SQMZ)*SGZANG   
        ASZIM=0.    
        TGZANG=1./SQMW*(BE2*(4.-2.*BE2+BE4)+BE2*(4.-10.*BE2+BE4)*CTH+   
     &  (2.-11.*BE2+10.*BE4)*CTH2+BE2*CTH3) 
        ATGRE=0.5*XW*SH/TH*TGZANG   
        ATGIM=0.    
        ATZRE=0.5*(1.-XW)*SH/(TH-SQMZ)*TGZANG   
        ATZIM=0.    
        A4RE=1./SQMW*(1.+2.*BE2-6.*BE2*CTH-CTH2)    
        A4IM=0. 
        FACH=COMFAC*1./(4096.*PARU(1)**2)*(AEM/XW)**4*(SH/SQMW)**2* 
     &  ((ASHRE+ATHRE+ASGRE+ASZRE+ATGRE+ATZRE+A4RE)**2+ 
     &  (ASHIM+ATHIM+ASGIM+ASZIM+ATGIM+ATZIM+A4IM)**2)  
        DO 610 I=MIN1,MAX1  
        IF(I.EQ.0.OR.KFAC(1,I).EQ.0) GOTO 610   
        EI=SIGN(1.,FLOAT(I))*KCHG(IABS(I),1)    
        DO 600 J=MIN2,MAX2  
        IF(J.EQ.0.OR.KFAC(2,J).EQ.0) GOTO 600   
        EJ=SIGN(1.,FLOAT(J))*KCHG(IABS(J),1)    
        IF(EI*EJ.GT.0.) GOTO 600    
        NCHN=NCHN+1 
        ISIG(NCHN,1)=I  
        ISIG(NCHN,2)=J  
        ISIG(NCHN,3)=1  
        SIGH(NCHN)=FACH*VINT(180+I)*VINT(180+J) 
  600   CONTINUE    
  610   CONTINUE    
    
      ELSEIF(ISUB.EQ.78) THEN   
C...W+/- + H0 -> W+/- + H0. 
    
      ELSEIF(ISUB.EQ.79) THEN   
C...H0 + H0 -> H0 + H0. 
    
      ENDIF 
    
C...C: 2 -> 2, tree diagrams with masses.   
    
      ELSEIF(ISUB.LE.90) THEN   
      IF(ISUB.EQ.81) THEN   
C...q + qb -> Q + QB.   
        FACQQB=COMFAC*AS**2*4./9.*(((TH-SQM3)**2+   
     &  (UH-SQM3)**2)/SH2+2.*SQM3/SH)   
        IF(MSTP(35).GE.1) THEN  
          IF(MSTP(35).EQ.1) THEN    
            ALSSG=PARP(35)  
          ELSE  
            MST115=MSTU(115)    
            MSTU(115)=MSTP(36)  
            Q2BN=SQRT(SQM3*((SQRT(SH)-2.*SQRT(SQM3))**2+PARP(36)**2))   
            ALSSG=ULALPS_HIJING(Q2BN)  
            MSTU(115)=MST115    
          ENDIF 
          XREPU=PARU(1)*ALSSG/(6.*SQRT(MAX(1E-20,1.-4.*SQM3/SH)))   
          FREPU=XREPU/(EXP(MIN(100.,XREPU))-1.) 
          PARI(81)=FREPU    
          FACQQB=FACQQB*FREPU   
        ENDIF   
        DO 620 I=MINA,MAXA  
        IF(I.EQ.0.OR.KFAC(1,I)*KFAC(2,-I).EQ.0) GOTO 620    
        NCHN=NCHN+1 
        ISIG(NCHN,1)=I  
        ISIG(NCHN,2)=-I 
        ISIG(NCHN,3)=1  
        SIGH(NCHN)=FACQQB   
  620   CONTINUE    
    
      ELSEIF(ISUB.EQ.82) THEN   
C...g + g -> Q + QB.    
        FACQQ1=COMFAC*FACA*AS**2*1./6.*((UH-SQM3)/(TH-SQM3)-    
     &  2.*(UH-SQM3)**2/SH2+4.*SQM3/SH*(TH*UH-SQM3**2)/(TH-SQM3)**2)    
        FACQQ2=COMFAC*FACA*AS**2*1./6.*((TH-SQM3)/(UH-SQM3)-    
     &  2.*(TH-SQM3)**2/SH2+4.*SQM3/SH*(TH*UH-SQM3**2)/(UH-SQM3)**2)    
        IF(MSTP(35).GE.1) THEN  
          IF(MSTP(35).EQ.1) THEN    
            ALSSG=PARP(35)  
          ELSE  
            MST115=MSTU(115)    
            MSTU(115)=MSTP(36)  
            Q2BN=SQRT(SQM3*((SQRT(SH)-2.*SQRT(SQM3))**2+PARP(36)**2))   
            ALSSG=ULALPS_HIJING(Q2BN)  
            MSTU(115)=MST115    
          ENDIF 
          XATTR=4.*PARU(1)*ALSSG/(3.*SQRT(MAX(1E-20,1.-4.*SQM3/SH)))    
          FATTR=XATTR/(1.-EXP(-MIN(100.,XATTR)))    
          XREPU=PARU(1)*ALSSG/(6.*SQRT(MAX(1E-20,1.-4.*SQM3/SH)))   
          FREPU=XREPU/(EXP(MIN(100.,XREPU))-1.) 
          FATRE=(2.*FATTR+5.*FREPU)/7.  
          PARI(81)=FATRE    
          FACQQ1=FACQQ1*FATRE   
          FACQQ2=FACQQ2*FATRE   
        ENDIF   
        IF(KFAC(1,21)*KFAC(2,21).EQ.0) GOTO 630 
        NCHN=NCHN+1 
        ISIG(NCHN,1)=21 
        ISIG(NCHN,2)=21 
        ISIG(NCHN,3)=1  
        SIGH(NCHN)=FACQQ1   
        NCHN=NCHN+1 
        ISIG(NCHN,1)=21 
        ISIG(NCHN,2)=21 
        ISIG(NCHN,3)=2  
        SIGH(NCHN)=FACQQ2   
  630   CONTINUE    
    
      ENDIF 
    
C...D: Mimimum bias processes.  
    
      ELSEIF(ISUB.LE.100) THEN  
      IF(ISUB.EQ.91) THEN   
C...Elastic scattering. 
        SIGS=XSEC(ISUB,1)   
    
      ELSEIF(ISUB.EQ.92) THEN   
C...Single diffractive scattering.  
        SIGS=XSEC(ISUB,1)   
    
      ELSEIF(ISUB.EQ.93) THEN   
C...Double diffractive scattering.  
        SIGS=XSEC(ISUB,1)   
    
      ELSEIF(ISUB.EQ.94) THEN   
C...Central diffractive scattering. 
        SIGS=XSEC(ISUB,1)   
    
      ELSEIF(ISUB.EQ.95) THEN   
C...Low-pT scattering.  
        SIGS=XSEC(ISUB,1)   
    
      ELSEIF(ISUB.EQ.96) THEN   
C...Multiple interactions: sum of QCD processes.    
        CALL PYWIDT_HIJING(21,SQRT(SH),WDTP,WDTE)  
    
C...q + q' -> q + q'.   
        FACQQ1=COMFAC*AS**2*4./9.*(SH2+UH2)/TH2 
        FACQQB=COMFAC*AS**2*4./9.*((SH2+UH2)/TH2*FACA-  
     &  MSTP(34)*2./3.*UH2/(SH*TH)) 
        FACQQ2=COMFAC*AS**2*4./9.*((SH2+TH2)/UH2-   
     &  MSTP(34)*2./3.*SH2/(TH*UH)) 
        DO 650 I=-3,3   
        IF(I.EQ.0) GOTO 650 
        DO 640 J=-3,3   
        IF(J.EQ.0) GOTO 640 
        NCHN=NCHN+1 
        ISIG(NCHN,1)=I  
        ISIG(NCHN,2)=J  
        ISIG(NCHN,3)=111    
        SIGH(NCHN)=FACQQ1   
        IF(I.EQ.-J) SIGH(NCHN)=FACQQB   
        IF(I.EQ.J) THEN 
          SIGH(NCHN)=0.5*SIGH(NCHN) 
          NCHN=NCHN+1   
          ISIG(NCHN,1)=I    
          ISIG(NCHN,2)=J    
          ISIG(NCHN,3)=112  
          SIGH(NCHN)=0.5*FACQQ2 
        ENDIF   
  640   CONTINUE    
  650   CONTINUE    
    
C...q + qb -> q' + qb' or g + g.    
        FACQQB=COMFAC*AS**2*4./9.*(TH2+UH2)/SH2*(WDTE(0,1)+WDTE(0,2)+   
     &  WDTE(0,3)+WDTE(0,4))    
        FACGG1=COMFAC*AS**2*32./27.*(UH/TH-(2.+MSTP(34)*1./4.)*UH2/SH2) 
        FACGG2=COMFAC*AS**2*32./27.*(TH/UH-(2.+MSTP(34)*1./4.)*TH2/SH2) 
        DO 660 I=-3,3   
        IF(I.EQ.0) GOTO 660 
        NCHN=NCHN+1 
        ISIG(NCHN,1)=I  
        ISIG(NCHN,2)=-I 
        ISIG(NCHN,3)=121    
        SIGH(NCHN)=FACQQB   
        NCHN=NCHN+1 
        ISIG(NCHN,1)=I  
        ISIG(NCHN,2)=-I 
        ISIG(NCHN,3)=131    
        SIGH(NCHN)=0.5*FACGG1   
        NCHN=NCHN+1 
        ISIG(NCHN,1)=I  
        ISIG(NCHN,2)=-I 
        ISIG(NCHN,3)=132    
        SIGH(NCHN)=0.5*FACGG2   
  660   CONTINUE    
    
C...q + g -> q + g. 
        FACQG1=COMFAC*AS**2*4./9.*((2.+MSTP(34)*1./4.)*UH2/TH2-UH/SH)*  
     &  FACA    
        FACQG2=COMFAC*AS**2*4./9.*((2.+MSTP(34)*1./4.)*SH2/TH2-SH/UH)   
        DO 680 I=-3,3   
        IF(I.EQ.0) GOTO 680 
        DO 670 ISDE=1,2 
        NCHN=NCHN+1 
        ISIG(NCHN,ISDE)=I   
        ISIG(NCHN,3-ISDE)=21    
        ISIG(NCHN,3)=281    
        SIGH(NCHN)=FACQG1   
        NCHN=NCHN+1 
        ISIG(NCHN,ISDE)=I   
        ISIG(NCHN,3-ISDE)=21    
        ISIG(NCHN,3)=282    
        SIGH(NCHN)=FACQG2   
  670   CONTINUE    
  680   CONTINUE    
    
C...g + g -> q + qb or g + g.   
        FACQQ1=COMFAC*AS**2*1./6.*(UH/TH-(2.+MSTP(34)*1./4.)*UH2/SH2)*  
     &  (WDTE(0,1)+WDTE(0,2)+WDTE(0,3)+WDTE(0,4))*FACA  
        FACQQ2=COMFAC*AS**2*1./6.*(TH/UH-(2.+MSTP(34)*1./4.)*TH2/SH2)*  
     &  (WDTE(0,1)+WDTE(0,2)+WDTE(0,3)+WDTE(0,4))*FACA  
        FACGG1=COMFAC*AS**2*9./4.*(SH2/TH2+2.*SH/TH+3.+2.*TH/SH+    
     &  TH2/SH2)*FACA   
        FACGG2=COMFAC*AS**2*9./4.*(UH2/SH2+2.*UH/SH+3.+2.*SH/UH+    
     &  SH2/UH2)*FACA   
        FACGG3=COMFAC*AS**2*9./4.*(TH2/UH2+2.*TH/UH+3+2.*UH/TH+UH2/TH2) 
        NCHN=NCHN+1 
        ISIG(NCHN,1)=21 
        ISIG(NCHN,2)=21 
        ISIG(NCHN,3)=531    
        SIGH(NCHN)=FACQQ1   
        NCHN=NCHN+1 
        ISIG(NCHN,1)=21 
        ISIG(NCHN,2)=21 
        ISIG(NCHN,3)=532    
        SIGH(NCHN)=FACQQ2   
        NCHN=NCHN+1 
        ISIG(NCHN,1)=21 
        ISIG(NCHN,2)=21 
        ISIG(NCHN,3)=681    
        SIGH(NCHN)=0.5*FACGG1   
        NCHN=NCHN+1 
        ISIG(NCHN,1)=21 
        ISIG(NCHN,2)=21 
        ISIG(NCHN,3)=682    
        SIGH(NCHN)=0.5*FACGG2   
        NCHN=NCHN+1 
        ISIG(NCHN,1)=21 
        ISIG(NCHN,2)=21 
        ISIG(NCHN,3)=683    
        SIGH(NCHN)=0.5*FACGG3   
      ENDIF 
    
C...E: 2 -> 1, loop diagrams.   
    
      ELSEIF(ISUB.LE.110) THEN  
      IF(ISUB.EQ.101) THEN  
C...g + g -> gamma*/Z0. 
    
      ELSEIF(ISUB.EQ.102) THEN  
C...g + g -> H0.    
        CALL PYWIDT_HIJING(25,SQRT(SH),WDTP,WDTE)  
        ETARE=0.    
        ETAIM=0.    
        DO 690 I=1,2*MSTP(1)    
        EPS=4.*PMAS(I,1)**2/SH  
        IF(EPS.LE.1.) THEN  
          IF(EPS.GT.1.E-4) THEN 
            ROOT=SQRT(1.-EPS)   
            RLN=LOG((1.+ROOT)/(1.-ROOT))    
          ELSE  
            RLN=LOG(4./EPS-2.)  
          ENDIF 
          PHIRE=0.25*(RLN**2-PARU(1)**2)    
          PHIIM=0.5*PARU(1)*RLN 
        ELSE    
          PHIRE=-(ASIN(1./SQRT(EPS)))**2    
          PHIIM=0.  
        ENDIF   
        ETARE=ETARE+0.5*EPS*(1.+(EPS-1.)*PHIRE) 
        ETAIM=ETAIM+0.5*EPS*(EPS-1.)*PHIIM  
  690   CONTINUE    
        ETA2=ETARE**2+ETAIM**2  
        FACH=COMFAC*FACA*(AS/PARU(1)*AEM/XW)**2*1./512.*    
     &  (SH/SQMW)**2*ETA2*SH2/((SH-SQMH)**2+GMMH**2)*   
     &  (WDTE(0,1)+WDTE(0,2)+WDTE(0,4)) 
        IF(KFAC(1,21)*KFAC(2,21).EQ.0) GOTO 700 
        NCHN=NCHN+1 
        ISIG(NCHN,1)=21 
        ISIG(NCHN,2)=21 
        ISIG(NCHN,3)=1  
        SIGH(NCHN)=FACH 
  700   CONTINUE    
    
      ENDIF 
    
C...F: 2 -> 2, box diagrams.    
    
      ELSEIF(ISUB.LE.120) THEN  
      IF(ISUB.EQ.111) THEN  
C...f + fb -> g + H0 (q + qb -> g + H0 only).   
        A5STUR=0.   
        A5STUI=0.   
        DO 710 I=1,2*MSTP(1)    
        SQMQ=PMAS(I,1)**2   
        EPSS=4.*SQMQ/SH 
        EPSH=4.*SQMQ/SQMH   
        A5STUR=A5STUR+SQMQ/SQMH*(4.+4.*SH/(TH+UH)*(PYW1AU_HIJING(EPSS,1)
     $       -PYW1AU_HIJING(EPSH,1))+(1.-4.*SQMQ/(TH+UH))
     $       *(PYW2AU_HIJING(EPSS,1)-PYW2AU_HIJING(EPSH,1)))    
        A5STUI=A5STUI+SQMQ/SQMH*(4.*SH/(TH+UH)*(PYW1AU_HIJING(EPSS,2)- 
     &       PYW1AU_HIJING(EPSH,2))+(1.-4.*SQMQ/(TH+UH))
     $       *(PYW2AU_HIJING(EPSS,2)-PYW2AU_HIJING(EPSH,2)))    
  710   CONTINUE    
        FACGH=COMFAC*FACA/(144.*PARU(1)**2)*AEM/XW*AS**3*SQMH/SQMW* 
     &  SQMH/SH*(UH**2+TH**2)/(UH+TH)**2*(A5STUR**2+A5STUI**2)  
        FACGH=FACGH*WIDS(25,2)  
        DO 720 I=MINA,MAXA  
        IF(I.EQ.0.OR.KFAC(1,I)*KFAC(2,-I).EQ.0) GOTO 720    
        NCHN=NCHN+1 
        ISIG(NCHN,1)=I  
        ISIG(NCHN,2)=-I 
        ISIG(NCHN,3)=1  
        SIGH(NCHN)=FACGH    
  720   CONTINUE    
    
      ELSEIF(ISUB.EQ.112) THEN  
C...f + g -> f + H0 (q + g -> q + H0 only). 
        A5TSUR=0.   
        A5TSUI=0.   
        DO 730 I=1,2*MSTP(1)    
        SQMQ=PMAS(I,1)**2   
        EPST=4.*SQMQ/TH 
        EPSH=4.*SQMQ/SQMH   
        A5TSUR=A5TSUR+SQMQ/SQMH*(4.+4.*TH/(SH+UH)*(PYW1AU_HIJING(EPST,1)
     $       -PYW1AU_HIJING(EPSH,1))+(1.-4.*SQMQ/(SH+UH))
     $       *(PYW2AU_HIJING(EPST,1)-PYW2AU_HIJING(EPSH,1)))    
        A5TSUI=A5TSUI+SQMQ/SQMH*(4.*TH/(SH+UH)*(PYW1AU_HIJING(EPST,2)- 
     &       PYW1AU_HIJING(EPSH,2))+(1.-4.*SQMQ/(SH+UH))
     $       *(PYW2AU_HIJING(EPST,2)-PYW2AU_HIJING(EPSH,2)))    
  730   CONTINUE    
        FACQH=COMFAC*FACA/(384.*PARU(1)**2)*AEM/XW*AS**3*SQMH/SQMW* 
     &  SQMH/(-TH)*(UH**2+SH**2)/(UH+SH)**2*(A5TSUR**2+A5TSUI**2)   
        FACQH=FACQH*WIDS(25,2)  
        DO 750 I=MINA,MAXA  
        IF(I.EQ.0) GOTO 750 
        DO 740 ISDE=1,2 
        IF(ISDE.EQ.1.AND.KFAC(1,I)*KFAC(2,21).EQ.0) GOTO 740    
        IF(ISDE.EQ.2.AND.KFAC(1,21)*KFAC(2,I).EQ.0) GOTO 740    
        NCHN=NCHN+1 
        ISIG(NCHN,ISDE)=I   
        ISIG(NCHN,3-ISDE)=21    
        ISIG(NCHN,3)=1  
        SIGH(NCHN)=FACQH    
  740   CONTINUE    
  750   CONTINUE    
    
      ELSEIF(ISUB.EQ.113) THEN  
C...g + g -> g + H0.    
        A2STUR=0.   
        A2STUI=0.   
        A2USTR=0.   
        A2USTI=0.   
        A2TUSR=0.   
        A2TUSI=0.   
        A4STUR=0.   
        A4STUI=0.   
        DO 760 I=6,2*MSTP(1)    
C'''Only t-quarks yet included  
        SQMQ=PMAS(I,1)**2   
        EPSS=4.*SQMQ/SH 
        EPST=4.*SQMQ/TH 
        EPSU=4.*SQMQ/UH 
        EPSH=4.*SQMQ/SQMH   
        IF(EPSH.LT.1.E-6) GOTO 760  
        BESTU=0.5*(1.+SQRT(1.+EPSS*TH/UH))  
        BEUST=0.5*(1.+SQRT(1.+EPSU*SH/TH))  
        BETUS=0.5*(1.+SQRT(1.+EPST*UH/SH))  
        BEUTS=BESTU 
        BETSU=BEUST 
        BESUT=BETUS 
        W3STUR=PYI3AU_HIJING(BESTU,EPSH,1)-PYI3AU_HIJING(BESTU,EPSS,1)-   
     &  PYI3AU_HIJING(BESTU,EPSU,1)    
        W3STUI=PYI3AU_HIJING(BESTU,EPSH,2)-PYI3AU_HIJING(BESTU,EPSS,2)-   
     &  PYI3AU_HIJING(BESTU,EPSU,2)    
        W3SUTR=PYI3AU_HIJING(BESUT,EPSH,1)-PYI3AU_HIJING(BESUT,EPSS,1)-   
     &  PYI3AU_HIJING(BESUT,EPST,1)    
        W3SUTI=PYI3AU_HIJING(BESUT,EPSH,2)-PYI3AU_HIJING(BESUT,EPSS,2)-   
     &  PYI3AU_HIJING(BESUT,EPST,2)    
        W3TSUR=PYI3AU_HIJING(BETSU,EPSH,1)-PYI3AU_HIJING(BETSU,EPST,1)-   
     &  PYI3AU_HIJING(BETSU,EPSU,1)    
        W3TSUI=PYI3AU_HIJING(BETSU,EPSH,2)-PYI3AU_HIJING(BETSU,EPST,2)-   
     &  PYI3AU_HIJING(BETSU,EPSU,2)    
        W3TUSR=PYI3AU_HIJING(BETUS,EPSH,1)-PYI3AU_HIJING(BETUS,EPST,1)-   
     &  PYI3AU_HIJING(BETUS,EPSS,1)    
        W3TUSI=PYI3AU_HIJING(BETUS,EPSH,2)-PYI3AU_HIJING(BETUS,EPST,2)-   
     &  PYI3AU_HIJING(BETUS,EPSS,2)    
        W3USTR=PYI3AU_HIJING(BEUST,EPSH,1)-PYI3AU_HIJING(BEUST,EPSU,1)-   
     &  PYI3AU_HIJING(BEUST,EPST,1)    
        W3USTI=PYI3AU_HIJING(BEUST,EPSH,2)-PYI3AU_HIJING(BEUST,EPSU,2)-   
     &  PYI3AU_HIJING(BEUST,EPST,2)    
        W3UTSR=PYI3AU_HIJING(BEUTS,EPSH,1)-PYI3AU_HIJING(BEUTS,EPSU,1)-   
     &  PYI3AU_HIJING(BEUTS,EPSS,1)    
        W3UTSI=PYI3AU_HIJING(BEUTS,EPSH,2)-PYI3AU_HIJING(BEUTS,EPSU,2)-   
     &  PYI3AU_HIJING(BEUTS,EPSS,2)    
        B2STUR=SQMQ/SQMH**2*(SH*(UH-SH)/(SH+UH)+2.*TH*UH*(UH+2.*SH)/    
     &       (SH+UH)**2*(PYW1AU_HIJING(EPST,1)-PYW1AU_HIJING(EPSH,1))
     $       +(SQMQ-SH/4.)*(0.5*PYW2AU_HIJING(EPSS,1)+0.5
     $       *PYW2AU_HIJING(EPSH,1)-PYW2AU_HIJING(EPST,1)+W3STUR)+SH**2
     $       *(2.*SQMQ/(SH+UH)**2-0.5/(SH+UH))*(PYW2AU_HIJING(EPST,1)
     $       -PYW2AU_HIJING(EPSH,1))+0.5*TH*UH/SH*(PYW2AU_HIJING(EPSH,1)
     $       -2.*PYW2AU_HIJING(EPST,1))+0.125*(SH-12.*SQMQ-4.*TH*UH/SH)
     $       *W3TSUR) 
        B2STUI=SQMQ/SQMH**2*(2.*TH*UH*(UH+2.*SH)/(SH+UH)**2*    
     &  (PYW1AU_HIJING(EPST,2)-PYW1AU_HIJING(EPSH,2))+(SQMQ-SH/4.)*   
     &       (0.5*PYW2AU_HIJING(EPSS,2)+0.5*PYW2AU_HIJING(EPSH,2)
     $       -PYW2AU_HIJING(EPST,2)+W3STUI)+SH**2*(2.*SQMQ/(SH+UH)**2-0
     $       .5/(SH+UH))*(PYW2AU_HIJING(EPST,2)-PYW2AU_HIJING(EPSH,2))+0
     $       .5*TH*UH/SH*(PYW2AU_HIJING(EPSH,2)-2.*PYW2AU_HIJING(EPST,2)
     $       )+0.125*(SH-12.*SQMQ-4.*TH*UH/SH)*W3TSUI) 
        B2SUTR=SQMQ/SQMH**2*(SH*(TH-SH)/(SH+TH)+2.*UH*TH*(TH+2.*SH)/    
     &       (SH+TH)**2*(PYW1AU_HIJING(EPSU,1)-PYW1AU_HIJING(EPSH,1))
     $       +(SQMQ-SH/4.)*(0.5*PYW2AU_HIJING(EPSS,1)+0.5
     $       *PYW2AU_HIJING(EPSH,1)-PYW2AU_HIJING(EPSU,1)+W3SUTR)+SH**2
     $       *(2.*SQMQ/(SH+TH)**2-0.5/(SH+TH))*(PYW2AU_HIJING(EPSU,1)
     $       -PYW2AU_HIJING(EPSH,1))+0.5*UH*TH/SH*(PYW2AU_HIJING(EPSH,1)
     $       -2.*PYW2AU_HIJING(EPSU,1))+0.125*(SH-12.*SQMQ-4.*UH*TH/SH)
     $       *W3USTR) 
        B2SUTI=SQMQ/SQMH**2*(2.*UH*TH*(TH+2.*SH)/(SH+TH)**2*    
     &  (PYW1AU_HIJING(EPSU,2)-PYW1AU_HIJING(EPSH,2))+(SQMQ-SH/4.)*   
     &       (0.5*PYW2AU_HIJING(EPSS,2)+0.5*PYW2AU_HIJING(EPSH,2)
     $       -PYW2AU_HIJING(EPSU,2)+W3SUTI)+SH**2*(2.*SQMQ/(SH+TH)**2-0
     $       .5/(SH+TH))*(PYW2AU_HIJING(EPSU,2)-PYW2AU_HIJING(EPSH,2))+0
     $       .5*UH*TH/SH*(PYW2AU_HIJING(EPSH,2)-2.*PYW2AU_HIJING(EPSU,2)
     $       )+0.125*(SH-12.*SQMQ-4.*UH*TH/SH)*W3USTI) 
        B2TSUR=SQMQ/SQMH**2*(TH*(UH-TH)/(TH+UH)+2.*SH*UH*(UH+2.*TH)/    
     &       (TH+UH)**2*(PYW1AU_HIJING(EPSS,1)-PYW1AU_HIJING(EPSH,1))
     $       +(SQMQ-TH/4.)*(0.5*PYW2AU_HIJING(EPST,1)+0.5
     $       *PYW2AU_HIJING(EPSH,1)-PYW2AU_HIJING(EPSS,1)+W3TSUR)+TH**2
     $       *(2.*SQMQ/(TH+UH)**2-0.5/(TH+UH))*(PYW2AU_HIJING(EPSS,1)
     $       -PYW2AU_HIJING(EPSH,1))+0.5*SH*UH/TH*(PYW2AU_HIJING(EPSH,1)
     $       -2.*PYW2AU_HIJING(EPSS,1))+0.125*(TH-12.*SQMQ-4.*SH*UH/TH)
     $       *W3STUR) 
        B2TSUI=SQMQ/SQMH**2*(2.*SH*UH*(UH+2.*TH)/(TH+UH)**2*    
     &  (PYW1AU_HIJING(EPSS,2)-PYW1AU_HIJING(EPSH,2))+(SQMQ-TH/4.)*   
     &       (0.5*PYW2AU_HIJING(EPST,2)+0.5*PYW2AU_HIJING(EPSH,2)
     $       -PYW2AU_HIJING(EPSS,2)+W3TSUI)+TH**2*(2.*SQMQ/(TH+UH)**2-0
     $       .5/(TH+UH))*(PYW2AU_HIJING(EPSS,2)-PYW2AU_HIJING(EPSH,2))+0
     $       .5*SH*UH/TH*(PYW2AU_HIJING(EPSH,2)-2.*PYW2AU_HIJING(EPSS,2)
     $       )+0.125*(TH-12.*SQMQ-4.*SH*UH/TH)*W3STUI) 
        B2TUSR=SQMQ/SQMH**2*(TH*(SH-TH)/(TH+SH)+2.*UH*SH*(SH+2.*TH)/    
     &       (TH+SH)**2*(PYW1AU_HIJING(EPSU,1)-PYW1AU_HIJING(EPSH,1))
     $       +(SQMQ-TH/4.)*(0.5*PYW2AU_HIJING(EPST,1)+0.5
     $       *PYW2AU_HIJING(EPSH,1)-PYW2AU_HIJING(EPSU,1)+W3TUSR)+TH**2
     $       *(2.*SQMQ/(TH+SH)**2-0.5/(TH+SH))*(PYW2AU_HIJING(EPSU,1)
     $       -PYW2AU_HIJING(EPSH,1))+0.5*UH*SH/TH*(PYW2AU_HIJING(EPSH,1)
     $       -2.*PYW2AU_HIJING(EPSU,1))+0.125*(TH-12.*SQMQ-4.*UH*SH/TH)
     $       *W3UTSR) 
        B2TUSI=SQMQ/SQMH**2*(2.*UH*SH*(SH+2.*TH)/(TH+SH)**2*    
     &  (PYW1AU_HIJING(EPSU,2)-PYW1AU_HIJING(EPSH,2))+(SQMQ-TH/4.)*   
     &       (0.5*PYW2AU_HIJING(EPST,2)+0.5*PYW2AU_HIJING(EPSH,2)
     $       -PYW2AU_HIJING(EPSU,2)+W3TUSI)+TH**2*(2.*SQMQ/(TH+SH)**2-0
     $       .5/(TH+SH))*(PYW2AU_HIJING(EPSU,2)-PYW2AU_HIJING(EPSH,2))+0
     $       .5*UH*SH/TH*(PYW2AU_HIJING(EPSH,2)-2.*PYW2AU_HIJING(EPSU,2)
     $       )+0.125*(TH-12.*SQMQ-4.*UH*SH/TH)*W3UTSI) 
        B2USTR=SQMQ/SQMH**2*(UH*(TH-UH)/(UH+TH)+2.*SH*TH*(TH+2.*UH)/    
     &       (UH+TH)**2*(PYW1AU_HIJING(EPSS,1)-PYW1AU_HIJING(EPSH,1))
     $       +(SQMQ-UH/4.)*(0.5*PYW2AU_HIJING(EPSU,1)+0.5
     $       *PYW2AU_HIJING(EPSH,1)-PYW2AU_HIJING(EPSS,1)+W3USTR)+UH**2
     $       *(2.*SQMQ/(UH+TH)**2-0.5/(UH+TH))*(PYW2AU_HIJING(EPSS,1)
     $       -PYW2AU_HIJING(EPSH,1))+0.5*SH*TH/UH*(PYW2AU_HIJING(EPSH,1)
     $       -2.*PYW2AU_HIJING(EPSS,1))+0.125*(UH-12.*SQMQ-4.*SH*TH/UH)
     $       *W3SUTR) 
        B2USTI=SQMQ/SQMH**2*(2.*SH*TH*(TH+2.*UH)/(UH+TH)**2*    
     &  (PYW1AU_HIJING(EPSS,2)-PYW1AU_HIJING(EPSH,2))+(SQMQ-UH/4.)*   
     &       (0.5*PYW2AU_HIJING(EPSU,2)+0.5*PYW2AU_HIJING(EPSH,2)
     $       -PYW2AU_HIJING(EPSS,2)+W3USTI)+UH**2*(2.*SQMQ/(UH+TH)**2-0
     $       .5/(UH+TH))*(PYW2AU_HIJING(EPSS,2)-PYW2AU_HIJING(EPSH,2))+0
     $       .5*SH*TH/UH*(PYW2AU_HIJING(EPSH,2)-2.*PYW2AU_HIJING(EPSS,2)
     $       )+0.125*(UH-12.*SQMQ-4.*SH*TH/UH)*W3SUTI) 
        B2UTSR=SQMQ/SQMH**2*(UH*(SH-UH)/(UH+SH)+2.*TH*SH*(SH+2.*UH)/    
     &       (UH+SH)**2*(PYW1AU_HIJING(EPST,1)-PYW1AU_HIJING(EPSH,1))
     $       +(SQMQ-UH/4.)*(0.5*PYW2AU_HIJING(EPSU,1)+0.5
     $       *PYW2AU_HIJING(EPSH,1)-PYW2AU_HIJING(EPST,1)+W3UTSR)+UH**2
     $       *(2.*SQMQ/(UH+SH)**2-0.5/(UH+SH))*(PYW2AU_HIJING(EPST,1)
     $       -PYW2AU_HIJING(EPSH,1))+0.5*TH*SH/UH*(PYW2AU_HIJING(EPSH,1)
     $       -2.*PYW2AU_HIJING(EPST,1))+0.125*(UH-12.*SQMQ-4.*TH*SH/UH)
     $       *W3TUSR) 
        B2UTSI=SQMQ/SQMH**2*(2.*TH*SH*(SH+2.*UH)/(UH+SH)**2*    
     &  (PYW1AU_HIJING(EPST,2)-PYW1AU_HIJING(EPSH,2))+(SQMQ-UH/4.)*   
     &       (0.5*PYW2AU_HIJING(EPSU,2)+0.5*PYW2AU_HIJING(EPSH,2)
     $       -PYW2AU_HIJING(EPST,2)+W3UTSI)+UH**2*(2.*SQMQ/(UH+SH)**2-0
     $       .5/(UH+SH))*(PYW2AU_HIJING(EPST,2)-PYW2AU_HIJING(EPSH,2))+0
     $       .5*TH*SH/UH*(PYW2AU_HIJING(EPSH,2)-2.*PYW2AU_HIJING(EPST,2)
     $       )+0.125*(UH-12.*SQMQ-4.*TH*SH/UH)*W3TUSI) 
        B4STUR=SQMQ/SQMH*(-2./3.+(SQMQ/SQMH-1./4.)*(PYW2AU_HIJING(EPSS,1
     $       )-PYW2AU_HIJING(EPSH,1)+W3STUR)) 
        B4STUI=SQMQ/SQMH*(SQMQ/SQMH-1./4.)*(PYW2AU_HIJING(EPSS,2)- 
     &  PYW2AU_HIJING(EPSH,2)+W3STUI)  
        B4TUSR=SQMQ/SQMH*(-2./3.+(SQMQ/SQMH-1./4.)*(PYW2AU_HIJING(EPST,1
     $       )-PYW2AU_HIJING(EPSH,1)+W3TUSR)) 
        B4TUSI=SQMQ/SQMH*(SQMQ/SQMH-1./4.)*(PYW2AU_HIJING(EPST,2)- 
     &  PYW2AU_HIJING(EPSH,2)+W3TUSI)  
        B4USTR=SQMQ/SQMH*(-2./3.+(SQMQ/SQMH-1./4.)*(PYW2AU_HIJING(EPSU,1
     $       )-PYW2AU_HIJING(EPSH,1)+W3USTR)) 
        B4USTI=SQMQ/SQMH*(SQMQ/SQMH-1./4.)*(PYW2AU_HIJING(EPSU,2)- 
     &  PYW2AU_HIJING(EPSH,2)+W3USTI)  
        A2STUR=A2STUR+B2STUR+B2SUTR 
        A2STUI=A2STUI+B2STUI+B2SUTI 
        A2USTR=A2USTR+B2USTR+B2UTSR 
        A2USTI=A2USTI+B2USTI+B2UTSI 
        A2TUSR=A2TUSR+B2TUSR+B2TSUR 
        A2TUSI=A2TUSI+B2TUSI+B2TSUI 
        A4STUR=A4STUR+B4STUR+B4USTR+B4TUSR  
        A4STUI=A4STUI+B4STUI+B4USTI+B4TUSI  
  760   CONTINUE    
        FACGH=COMFAC*FACA*3./(128.*PARU(1)**2)*AEM/XW*AS**3*    
     &  SQMH/SQMW*SQMH**3/(SH*TH*UH)*(A2STUR**2+A2STUI**2+A2USTR**2+    
     &  A2USTI**2+A2TUSR**2+A2TUSI**2+A4STUR**2+A4STUI**2)  
        FACGH=FACGH*WIDS(25,2)  
        IF(KFAC(1,21)*KFAC(2,21).EQ.0) GOTO 770 
        NCHN=NCHN+1 
        ISIG(NCHN,1)=21 
        ISIG(NCHN,2)=21 
        ISIG(NCHN,3)=1  
        SIGH(NCHN)=FACGH    
  770   CONTINUE    
    
      ELSEIF(ISUB.EQ.114) THEN  
C...g + g -> gamma + gamma. 
        ASRE=0. 
        ASIM=0. 
        DO 780 I=1,2*MSTP(1)    
        EI=KCHG(IABS(I),1)/3.   
        SQMQ=PMAS(I,1)**2   
        EPSS=4.*SQMQ/SH 
        EPST=4.*SQMQ/TH 
        EPSU=4.*SQMQ/UH 
        IF(EPSS+ABS(EPST)+ABS(EPSU).LT.3.E-6) THEN  
          A0STUR=1.+(TH-UH)/SH*LOG(TH/UH)+0.5*(TH2+UH2)/SH2*    
     &    (LOG(TH/UH)**2+PARU(1)**2)    
          A0STUI=0. 
          A0TSUR=1.+(SH-UH)/TH*LOG(-SH/UH)+0.5*(SH2+UH2)/TH2*   
     &    LOG(-SH/UH)**2    
          A0TSUI=-PARU(1)*((SH-UH)/TH+(SH2+UH2)/TH2*LOG(-SH/UH))    
          A0UTSR=1.+(TH-SH)/UH*LOG(-TH/SH)+0.5*(TH2+SH2)/UH2*   
     &    LOG(-TH/SH)**2    
          A0UTSI=PARU(1)*((TH-SH)/UH+(TH2+SH2)/UH2*LOG(-TH/SH)) 
          A1STUR=-1.    
          A1STUI=0. 
          A2STUR=-1.    
          A2STUI=0. 
        ELSE    
          BESTU=0.5*(1.+SQRT(1.+EPSS*TH/UH))    
          BEUST=0.5*(1.+SQRT(1.+EPSU*SH/TH))    
          BETUS=0.5*(1.+SQRT(1.+EPST*UH/SH))    
          BEUTS=BESTU   
          BETSU=BEUST   
          BESUT=BETUS   
          A0STUR=1.+(1.+2.*TH/SH)*PYW1AU_HIJING(EPST,1)+(1.+2.*UH/SH)* 
     &         PYW1AU_HIJING(EPSU,1)+0.5*((TH2+UH2)/SH2-EPSS)
     $         *(PYW2AU_HIJING(EPST,1)+PYW2AU_HIJING(EPSU,1))-0.25*EPST
     $         *(1.-0.5*EPSS)*(PYI3AU_HIJING(BESUT,EPSS,1)
     $         +PYI3AU_HIJING(BESUT,EPST,1))-0.25*EPSU*(1.-0.5*EPSS)
     $         *(PYI3AU_HIJING(BESTU,EPSS,1)+PYI3AU_HIJING(BESTU,EPSU,1)
     $         )+0.25*(-2.*(TH2+UH2)/SH2+4.*EPSS+EPST+EPSU+0.5*EPST*EPSU
     $         )*(PYI3AU_HIJING(BETSU,EPST,1)+PYI3AU_HIJING(BETSU,EPSU,1
     $         ))   
          A0STUI=(1.+2.*TH/SH)*PYW1AU_HIJING(EPST,2)+(1.+2.*UH/SH)*    
     &         PYW1AU_HIJING(EPSU,2)+0.5*((TH2+UH2)/SH2-EPSS)
     $         *(PYW2AU_HIJING(EPST,2)+PYW2AU_HIJING(EPSU,2))-0.25*EPST
     $         *(1.-0.5*EPSS)*(PYI3AU_HIJING(BESUT,EPSS,2)
     $         +PYI3AU_HIJING(BESUT,EPST,2))-0.25*EPSU*(1.-0.5*EPSS)
     $         *(PYI3AU_HIJING(BESTU,EPSS,2)+PYI3AU_HIJING(BESTU,EPSU,2)
     $         )+0.25*(-2.*(TH2+UH2)/SH2+4.*EPSS+EPST+EPSU+0.5*EPST*EPSU
     $         )*(PYI3AU_HIJING(BETSU,EPST,2)+PYI3AU_HIJING(BETSU,EPSU,2
     $         ))   
          A0TSUR=1.+(1.+2.*SH/TH)*PYW1AU_HIJING(EPSS,1)+(1.+2.*UH/TH)* 
     &         PYW1AU_HIJING(EPSU,1)+0.5*((SH2+UH2)/TH2-EPST)
     $         *(PYW2AU_HIJING(EPSS,1)+PYW2AU_HIJING(EPSU,1))-0.25*EPSS
     $         *(1.-0.5*EPST)*(PYI3AU_HIJING(BETUS,EPST,1)
     $         +PYI3AU_HIJING(BETUS,EPSS,1))-0.25*EPSU*(1.-0.5*EPST)
     $         *(PYI3AU_HIJING(BETSU,EPST,1)+PYI3AU_HIJING(BETSU,EPSU,1)
     $         )+0.25*(-2.*(SH2+UH2)/TH2+4.*EPST+EPSS+EPSU+0.5*EPSS*EPSU
     $         )*(PYI3AU_HIJING(BESTU,EPSS,1)+PYI3AU_HIJING(BESTU,EPSU,1
     $         ))   
          A0TSUI=(1.+2.*SH/TH)*PYW1AU_HIJING(EPSS,2)+(1.+2.*UH/TH)*    
     &         PYW1AU_HIJING(EPSU,2)+0.5*((SH2+UH2)/TH2-EPST)
     $         *(PYW2AU_HIJING(EPSS,2)+PYW2AU_HIJING(EPSU,2))-0.25*EPSS
     $         *(1.-0.5*EPST)*(PYI3AU_HIJING(BETUS,EPST,2)
     $         +PYI3AU_HIJING(BETUS,EPSS,2))-0.25*EPSU*(1.-0.5*EPST)
     $         *(PYI3AU_HIJING(BETSU,EPST,2)+PYI3AU_HIJING(BETSU,EPSU,2)
     $         )+0.25*(-2.*(SH2+UH2)/TH2+4.*EPST+EPSS+EPSU+0.5*EPSS*EPSU
     $         )*(PYI3AU_HIJING(BESTU,EPSS,2)+PYI3AU_HIJING(BESTU,EPSU,2
     $         ))   
          A0UTSR=1.+(1.+2.*TH/UH)*PYW1AU_HIJING(EPST,1)+(1.+2.*SH/UH)* 
     &         PYW1AU_HIJING(EPSS,1)+0.5*((TH2+SH2)/UH2-EPSU)
     $         *(PYW2AU_HIJING(EPST,1)+PYW2AU_HIJING(EPSS,1))-0.25*EPST
     $         *(1.-0.5*EPSU)*(PYI3AU_HIJING(BEUST,EPSU,1)
     $         +PYI3AU_HIJING(BEUST,EPST,1))-0.25*EPSS*(1.-0.5*EPSU)
     $         *(PYI3AU_HIJING(BEUTS,EPSU,1)+PYI3AU_HIJING(BEUTS,EPSS,1)
     $         )+0.25*(-2.*(TH2+SH2)/UH2+4.*EPSU+EPST+EPSS+0.5*EPST*EPSS
     $         )*(PYI3AU_HIJING(BETUS,EPST,1)+PYI3AU_HIJING(BETUS,EPSS,1
     $         ))   
          A0UTSI=(1.+2.*TH/UH)*PYW1AU_HIJING(EPST,2)+(1.+2.*SH/UH)*    
     &         PYW1AU_HIJING(EPSS,2)+0.5*((TH2+SH2)/UH2-EPSU)
     $         *(PYW2AU_HIJING(EPST,2)+PYW2AU_HIJING(EPSS,2))-0.25*EPST
     $         *(1.-0.5*EPSU)*(PYI3AU_HIJING(BEUST,EPSU,2)
     $         +PYI3AU_HIJING(BEUST,EPST,2))-0.25*EPSS*(1.-0.5*EPSU)
     $         *(PYI3AU_HIJING(BEUTS,EPSU,2)+PYI3AU_HIJING(BEUTS,EPSS,2)
     $         )+0.25*(-2.*(TH2+SH2)/UH2+4.*EPSU+EPST+EPSS+0.5*EPST*EPSS
     $         )*(PYI3AU_HIJING(BETUS,EPST,2)+PYI3AU_HIJING(BETUS,EPSS,2
     $         ))   
          A1STUR=-1.-0.25*(EPSS+EPST+EPSU)*(PYW2AU_HIJING(EPSS,1)+ 
     &         PYW2AU_HIJING(EPST,1)+PYW2AU_HIJING(EPSU,1))+0.25*(EPSU+0
     $         .5*EPSS*EPST)*(PYI3AU_HIJING(BESUT,EPSS,1)
     $         +PYI3AU_HIJING(BESUT,EPST,1))+0.25*(EPST+0.5*EPSS*EPSU)
     $         *(PYI3AU_HIJING(BESTU,EPSS,1)+PYI3AU_HIJING(BESTU,EPSU,1)
     $         )+0.25*(EPSS+0.5*EPST*EPSU)*(PYI3AU_HIJING(BETSU,EPST,1)
     $         +PYI3AU_HIJING(BETSU,EPSU,1))   
          A1STUI=-0.25*(EPSS+EPST+EPSU)*(PYW2AU_HIJING(EPSS,2)
     $         +PYW2AU_HIJING(EPST,2)+PYW2AU_HIJING(EPSU,2))+0.25*(EPSU
     $         +0.5*EPSS*EPST)*(PYI3AU_HIJING(BESUT,EPSS,2)
     $         +PYI3AU_HIJING(BESUT,EPST,2))+0.25*(EPST+0.5*EPSS*EPSU)
     $         *(PYI3AU_HIJING(BESTU,EPSS,2)+PYI3AU_HIJING(BESTU,EPSU,2)
     $         )+0.25*(EPSS+0.5*EPST*EPSU)*(PYI3AU_HIJING(BETSU,EPST,2)
     $         +PYI3AU_HIJING(BETSU,EPSU,2))   
          A2STUR=-1.+0.125*EPSS*EPST*(PYI3AU_HIJING(BESUT,EPSS,1)+ 
     &         PYI3AU_HIJING(BESUT,EPST,1))+0.125*EPSS*EPSU
     $         *(PYI3AU_HIJING(BESTU,EPSS,1)+PYI3AU_HIJING(BESTU,EPSU,1)
     $         )+0.125*EPST*EPSU*(PYI3AU_HIJING(BETSU,EPST,1)
     $         +PYI3AU_HIJING(BETSU,EPSU,1)) 
          A2STUI=0.125*EPSS*EPST*(PYI3AU_HIJING(BESUT,EPSS,2)+ 
     &         PYI3AU_HIJING(BESUT,EPST,2))+0.125*EPSS*EPSU
     $         *(PYI3AU_HIJING(BESTU,EPSS,2)+PYI3AU_HIJING(BESTU,EPSU,2)
     $         )+0.125*EPST*EPSU*(PYI3AU_HIJING(BETSU,EPST,2)
     $         +PYI3AU_HIJING(BETSU,EPSU,2)) 
        ENDIF   
        ASRE=ASRE+EI**2*(A0STUR+A0TSUR+A0UTSR+4.*A1STUR+A2STUR) 
        ASIM=ASIM+EI**2*(A0STUI+A0TSUI+A0UTSI+4.*A1STUI+A2STUI) 
  780   CONTINUE    
        FACGG=COMFAC*FACA/(8.*PARU(1)**2)*AS**2*AEM**2*(ASRE**2+ASIM**2)    
        IF(KFAC(1,21)*KFAC(2,21).EQ.0) GOTO 790 
        NCHN=NCHN+1 
        ISIG(NCHN,1)=21 
        ISIG(NCHN,2)=21 
        ISIG(NCHN,3)=1  
        SIGH(NCHN)=FACGG    
  790   CONTINUE    
    
      ELSEIF(ISUB.EQ.115) THEN  
C...g + g -> gamma + Z0.    
    
      ELSEIF(ISUB.EQ.116) THEN  
C...g + g -> Z0 + Z0.   
    
      ELSEIF(ISUB.EQ.117) THEN  
C...g + g -> W+ + W-.   
    
      ENDIF 
    
C...G: 2 -> 3, tree diagrams.   
    
      ELSEIF(ISUB.LE.140) THEN  
      IF(ISUB.EQ.121) THEN  
C...g + g -> f + fb + H0.   
    
      ENDIF 
    
C...H: 2 -> 1, tree diagrams, non-standard model processes. 
    
      ELSEIF(ISUB.LE.160) THEN  
      IF(ISUB.EQ.141) THEN  
C...f + fb -> gamma*/Z0/Z'0.    
        MINT(61)=2  
        CALL PYWIDT_HIJING(32,SQRT(SH),WDTP,WDTE)  
        FACZP=COMFAC*AEM**2*4./9.   
        DO 800 I=MINA,MAXA  
        IF(I.EQ.0.OR.KFAC(1,I)*KFAC(2,-I).EQ.0) GOTO 800    
        EI=KCHG(IABS(I),1)/3.   
        AI=SIGN(1.,EI)  
        VI=AI-4.*EI*XW  
        API=SIGN(1.,EI) 
        VPI=API-4.*EI*XW    
        NCHN=NCHN+1 
        ISIG(NCHN,1)=I  
        ISIG(NCHN,2)=-I 
        ISIG(NCHN,3)=1  
        SIGH(NCHN)=FACZP*(EI**2*VINT(111)+EI*VI/(8.*XW*(1.-XW))*    
     &  SH*(SH-SQMZ)/((SH-SQMZ)**2+GMMZ**2)*VINT(112)+EI*VPI/(8.*XW*    
     &  (1.-XW))*SH*(SH-SQMZP)/((SH-SQMZP)**2+GMMZP**2)*VINT(113)+  
     &  (VI**2+AI**2)/(16.*XW*(1.-XW))**2*SH2/((SH-SQMZ)**2+GMMZ**2)*   
     &  VINT(114)+2.*(VI*VPI+AI*API)/(16.*XW*(1.-XW))**2*SH2*   
     &  ((SH-SQMZ)*(SH-SQMZP)+GMMZ*GMMZP)/(((SH-SQMZ)**2+GMMZ**2)*  
     &  ((SH-SQMZP)**2+GMMZP**2))*VINT(115)+(VPI**2+API**2)/    
     &  (16.*XW*(1.-XW))**2*SH2/((SH-SQMZP)**2+GMMZP**2)*VINT(116)) 
  800   CONTINUE    
    
      ELSEIF(ISUB.EQ.142) THEN  
C...f + fb' -> H+/-.    
        CALL PYWIDT_HIJING(37,SQRT(SH),WDTP,WDTE)  
        FHC=COMFAC*(AEM/XW)**2*1./48.*(SH/SQMW)**2*SH2/ 
     &  ((SH-SQMHC)**2+GMMHC**2)    
C'''No construction yet for leptons 
        DO 840 I=1,MSTP(54)/2   
        IL=2*I-1    
        IU=2*I  
        RMQL=PMAS(IL,1)**2/SH   
        RMQU=PMAS(IU,1)**2/SH   
        FACHC=FHC*((RMQL*PARU(121)+RMQU/PARU(121))*(1.-RMQL-RMQU)-  
     &  4.*RMQL*RMQU)/SQRT(MAX(0.,(1.-RMQL-RMQU)**2-4.*RMQL*RMQU))  
        IF(KFAC(1,IL)*KFAC(2,-IU).EQ.0) GOTO 810    
        KCHHC=(KCHG(IL,1)-KCHG(IU,1))/3 
        NCHN=NCHN+1 
        ISIG(NCHN,1)=IL 
        ISIG(NCHN,2)=-IU    
        ISIG(NCHN,3)=1  
        SIGH(NCHN)=FACHC*(WDTE(0,1)+WDTE(0,(5-KCHHC)/2)+WDTE(0,4))  
  810   IF(KFAC(1,-IL)*KFAC(2,IU).EQ.0) GOTO 820    
        KCHHC=(-KCHG(IL,1)+KCHG(IU,1))/3    
        NCHN=NCHN+1 
        ISIG(NCHN,1)=-IL    
        ISIG(NCHN,2)=IU 
        ISIG(NCHN,3)=1  
        SIGH(NCHN)=FACHC*(WDTE(0,1)+WDTE(0,(5-KCHHC)/2)+WDTE(0,4))  
  820   IF(KFAC(1,IU)*KFAC(2,-IL).EQ.0) GOTO 830    
        KCHHC=(KCHG(IU,1)-KCHG(IL,1))/3 
        NCHN=NCHN+1 
        ISIG(NCHN,1)=IU 
        ISIG(NCHN,2)=-IL    
        ISIG(NCHN,3)=1  
        SIGH(NCHN)=FACHC*(WDTE(0,1)+WDTE(0,(5-KCHHC)/2)+WDTE(0,4))  
  830   IF(KFAC(1,-IU)*KFAC(2,IL).EQ.0) GOTO 840    
        KCHHC=(-KCHG(IU,1)+KCHG(IL,1))/3    
        NCHN=NCHN+1 
        ISIG(NCHN,1)=-IU    
        ISIG(NCHN,2)=IL 
        ISIG(NCHN,3)=1  
        SIGH(NCHN)=FACHC*(WDTE(0,1)+WDTE(0,(5-KCHHC)/2)+WDTE(0,4))  
  840   CONTINUE    
    
      ELSEIF(ISUB.EQ.143) THEN  
C...f + fb -> R.    
        CALL PYWIDT_HIJING(40,SQRT(SH),WDTP,WDTE)  
        FACR=COMFAC*(AEM/XW)**2*1./9.*SH2/((SH-SQMR)**2+GMMR**2)    
        DO 860 I=MIN1,MAX1  
        IF(I.EQ.0.OR.KFAC(1,I).EQ.0) GOTO 860   
        IA=IABS(I)  
        DO 850 J=MIN2,MAX2  
        IF(J.EQ.0.OR.KFAC(2,J).EQ.0) GOTO 850   
        JA=IABS(J)  
        IF(I*J.GT.0.OR.IABS(IA-JA).NE.2) GOTO 850   
        NCHN=NCHN+1 
        ISIG(NCHN,1)=I  
        ISIG(NCHN,2)=J  
        ISIG(NCHN,3)=1  
        SIGH(NCHN)=FACR*(WDTE(0,1)+WDTE(0,(10-(I+J))/4)+WDTE(0,4))  
  850   CONTINUE    
  860   CONTINUE    
    
      ENDIF 
    
C...I: 2 -> 2, tree diagrams, non-standard model processes. 
    
      ELSE  
      IF(ISUB.EQ.161) THEN  
C...f + g -> f' + H+/- (q + g -> q' + H+/- only).   
        FHCQ=COMFAC*FACA*AS*AEM/XW*1./24    
        DO 900 I=1,MSTP(54) 
        IU=I+MOD(I,2)   
        SQMQ=PMAS(IU,1)**2  
        FACHCQ=FHCQ/PARU(121)*SQMQ/SQMW*(SH/(SQMQ-UH)+  
     &  2.*SQMQ*(SQMHC-UH)/(SQMQ-UH)**2+(SQMQ-UH)/SH+   
     &  2.*SQMQ/(SQMQ-UH)+2.*(SQMHC-UH)/(SQMQ-UH)*(SQMHC-SQMQ-SH)/SH)   
        IF(KFAC(1,-I)*KFAC(2,21).EQ.0) GOTO 870 
        KCHHC=ISIGN(1,-KCHG(I,1))   
        NCHN=NCHN+1 
        ISIG(NCHN,1)=-I 
        ISIG(NCHN,2)=21 
        ISIG(NCHN,3)=1  
        SIGH(NCHN)=FACHCQ*(WDTE(0,1)+WDTE(0,(5-KCHHC)/2)+WDTE(0,4)) 
  870   IF(KFAC(1,I)*KFAC(2,21).EQ.0) GOTO 880  
        KCHHC=ISIGN(1,KCHG(I,1))    
        NCHN=NCHN+1 
        ISIG(NCHN,1)=I  
        ISIG(NCHN,2)=21 
        ISIG(NCHN,3)=1  
        SIGH(NCHN)=FACHCQ*(WDTE(0,1)+WDTE(0,(5-KCHHC)/2)+WDTE(0,4)) 
  880   IF(KFAC(1,21)*KFAC(2,-I).EQ.0) GOTO 890 
        KCHHC=ISIGN(1,-KCHG(I,1))   
        NCHN=NCHN+1 
        ISIG(NCHN,1)=21 
        ISIG(NCHN,2)=-I 
        ISIG(NCHN,3)=1  
        SIGH(NCHN)=FACHCQ*(WDTE(0,1)+WDTE(0,(5-KCHHC)/2)+WDTE(0,4)) 
  890   IF(KFAC(1,21)*KFAC(2,I).EQ.0) GOTO 900  
        KCHHC=ISIGN(1,KCHG(I,1))    
        NCHN=NCHN+1 
        ISIG(NCHN,1)=21 
        ISIG(NCHN,2)=I  
        ISIG(NCHN,3)=1  
        SIGH(NCHN)=FACHCQ*(WDTE(0,1)+WDTE(0,(5-KCHHC)/2)+WDTE(0,4)) 
  900   CONTINUE    
    
      ENDIF 
      ENDIF 
    
C...Multiply with structure functions.  
      IF(ISUB.LE.90.OR.ISUB.GE.96) THEN 
        DO 910 ICHN=1,NCHN  
        IF(MINT(41).EQ.2) THEN  
          KFL1=ISIG(ICHN,1) 
          IF(KFL1.EQ.21) KFL1=0 
          SIGH(ICHN)=SIGH(ICHN)*XSFX(1,KFL1)    
        ENDIF   
        IF(MINT(42).EQ.2) THEN  
          KFL2=ISIG(ICHN,2) 
          IF(KFL2.EQ.21) KFL2=0 
          SIGH(ICHN)=SIGH(ICHN)*XSFX(2,KFL2)    
        ENDIF   
  910   SIGS=SIGS+SIGH(ICHN)    
      ENDIF 
    
      RETURN    
      END   
