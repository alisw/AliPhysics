* $Id$
    
C*********************************************************************  
    
      SUBROUTINE PYRAND_HIJING 
    
C...Generates quantities characterizing the high-pT scattering at the   
C...parton level according to the matrix elements. Chooses incoming,    
C...reacting partons, their momentum fractions and one of the possible  
C...subprocesses.   
#include "ludat1_hijing.inc"
#include "ludat2_hijing.inc"
#include "pysubs_hijing.inc"
#include "pypars_hijing.inc"
#include "pyint1_hijing.inc"
#include "pyint2_hijing.inc"
#include "pyint3_hijing.inc"
#include "pyint4_hijing.inc"
#include "pyint5_hijing.inc"
    
C...Initial values, specifically for (first) semihard interaction.  
      MINT(17)=0    
      MINT(18)=0    
      VINT(143)=1.  
      VINT(144)=1.  
      IF(MSUB(95).EQ.1.OR.MINT(82).GE.2) CALL PYMULT_HIJING(2) 
      ISUB=0    
  100 MINT(51)=0    
    
C...Choice of process type - first event of overlay.    
      IF(MINT(82).EQ.1.AND.(ISUB.LE.90.OR.ISUB.GT.96)) THEN 
        RSUB=XSEC(0,1)*RLU_HIJING(0)   
        DO 110 I=1,200  
        IF(MSUB(I).NE.1) GOTO 110   
        ISUB=I  
        RSUB=RSUB-XSEC(I,1) 
        IF(RSUB.LE.0.) GOTO 120 
  110   CONTINUE    
  120   IF(ISUB.EQ.95) ISUB=96  
    
C...Choice of inclusive process type - overlayed events.    
      ELSEIF(MINT(82).GE.2.AND.ISUB.EQ.0) THEN  
        RSUB=VINT(131)*RLU_HIJING(0)   
        ISUB=96 
        IF(RSUB.GT.VINT(106)) ISUB=93   
        IF(RSUB.GT.VINT(106)+VINT(104)) ISUB=92 
        IF(RSUB.GT.VINT(106)+VINT(104)+VINT(103)) ISUB=91   
      ENDIF 
      IF(MINT(82).EQ.1) NGEN(0,1)=NGEN(0,1)+1   
      IF(MINT(82).EQ.1) NGEN(ISUB,1)=NGEN(ISUB,1)+1 
      MINT(1)=ISUB  
    
C...Find resonances (explicit or implicit in cross-section).    
      MINT(72)=0    
      KFR1=0    
      IF(ISET(ISUB).EQ.1.OR.ISET(ISUB).EQ.3) THEN   
        KFR1=KFPR(ISUB,1)   
      ELSEIF(ISUB.GE.71.AND.ISUB.LE.77) THEN    
        KFR1=25 
      ENDIF 
      IF(KFR1.NE.0) THEN    
        TAUR1=PMAS(KFR1,1)**2/VINT(2)   
        GAMR1=PMAS(KFR1,1)*PMAS(KFR1,2)/VINT(2) 
        MINT(72)=1  
        MINT(73)=KFR1   
        VINT(73)=TAUR1  
        VINT(74)=GAMR1  
      ENDIF 
      IF(ISUB.EQ.141) THEN  
        KFR2=23 
        TAUR2=PMAS(KFR2,1)**2/VINT(2)   
        GAMR2=PMAS(KFR2,1)*PMAS(KFR2,2)/VINT(2) 
        MINT(72)=2  
        MINT(74)=KFR2   
        VINT(75)=TAUR2  
        VINT(76)=GAMR2  
      ENDIF 
    
C...Find product masses and minimum pT of process,  
C...optionally with broadening according to a truncated Breit-Wigner.   
      VINT(63)=0.   
      VINT(64)=0.   
      MINT(71)=0    
      VINT(71)=CKIN(3)  
      IF(MINT(82).GE.2) VINT(71)=0. 
      IF(ISET(ISUB).EQ.2.OR.ISET(ISUB).EQ.4) THEN   
        DO 130 I=1,2    
        IF(KFPR(ISUB,I).EQ.0) THEN  
        ELSEIF(MSTP(42).LE.0) THEN  
          VINT(62+I)=PMAS(KFPR(ISUB,I),1)**2    
        ELSE    
          VINT(62+I)=ULMASS_HIJING(KFPR(ISUB,I))**2    
        ENDIF   
  130   CONTINUE    
        IF(MIN(VINT(63),VINT(64)).LT.CKIN(6)**2) MINT(71)=1 
        IF(MINT(71).EQ.1) VINT(71)=MAX(CKIN(3),CKIN(5)) 
      ENDIF 
    
      IF(ISET(ISUB).EQ.0) THEN  
C...Double or single diffractive, or elastic scattering:    
C...choose m^2 according to 1/m^2 (diffractive), constant (elastic) 
        IS=INT(1.5+RLU_HIJING(0))  
        VINT(63)=VINT(3)**2 
        VINT(64)=VINT(4)**2 
        IF(ISUB.EQ.92.OR.ISUB.EQ.93) VINT(62+IS)=PARP(111)**2   
        IF(ISUB.EQ.93) VINT(65-IS)=PARP(111)**2 
        SH=VINT(2)  
        SQM1=VINT(3)**2 
        SQM2=VINT(4)**2 
        SQM3=VINT(63)   
        SQM4=VINT(64)   
        SQLA12=(SH-SQM1-SQM2)**2-4.*SQM1*SQM2   
        SQLA34=(SH-SQM3-SQM4)**2-4.*SQM3*SQM4   
        THTER1=SQM1+SQM2+SQM3+SQM4-(SQM1-SQM2)*(SQM3-SQM4)/SH-SH    
        THTER2=SQRT(MAX(0.,SQLA12))*SQRT(MAX(0.,SQLA34))/SH 
        THL=0.5*(THTER1-THTER2) 
        THU=0.5*(THTER1+THTER2) 
        THM=MIN(MAX(THL,PARP(101)),THU) 
        JTMAX=0 
        IF(ISUB.EQ.92.OR.ISUB.EQ.93) JTMAX=ISUB-91  
        DO 140 JT=1,JTMAX   
        MINT(13+3*JT-IS*(2*JT-3))=1 
        SQMMIN=VINT(59+3*JT-IS*(2*JT-3))    
        SQMI=VINT(8-3*JT+IS*(2*JT-3))**2    
        SQMJ=VINT(3*JT-1-IS*(2*JT-3))**2    
        SQMF=VINT(68-3*JT+IS*(2*JT-3))  
        SQUA=0.5*SH/SQMI*((1.+(SQMI-SQMJ)/SH)*THM+SQMI-SQMF-    
     &  SQMJ**2/SH+(SQMI+SQMJ)*SQMF/SH+(SQMI-SQMJ)**2/SH**2*SQMF)   
        QUAR=SH/SQMI*(THM*(THM+SH-SQMI-SQMJ-SQMF*(1.-(SQMI-SQMJ)/SH))+  
     &  SQMI*SQMJ-SQMJ*SQMF*(1.+(SQMI-SQMJ-SQMF)/SH))   
        SQMMAX=SQUA+SQRT(MAX(0.,SQUA**2-QUAR))  
        IF(ABS(QUAR/SQUA**2).LT.1.E-06) SQMMAX=0.5*QUAR/SQUA    
        SQMMAX=MIN(SQMMAX,(VINT(1)-SQRT(SQMF))**2)  
        VINT(59+3*JT-IS*(2*JT-3))=SQMMIN*(SQMMAX/SQMMIN)**RLU_HIJING(0)    
  140   CONTINUE    
C...Choose t-hat according to exp(B*t-hat+C*t-hat^2).   
        SQM3=VINT(63)   
        SQM4=VINT(64)   
        SQLA34=(SH-SQM3-SQM4)**2-4.*SQM3*SQM4   
        THTER1=SQM1+SQM2+SQM3+SQM4-(SQM1-SQM2)*(SQM3-SQM4)/SH-SH    
        THTER2=SQRT(MAX(0.,SQLA12))*SQRT(MAX(0.,SQLA34))/SH 
        THL=0.5*(THTER1-THTER2) 
        THU=0.5*(THTER1+THTER2) 
        B=VINT(121) 
        C=VINT(122) 
        IF(ISUB.EQ.92.OR.ISUB.EQ.93) THEN   
          B=0.5*B   
          C=0.5*C   
        ENDIF   
        THM=MIN(MAX(THL,PARP(101)),THU) 
        EXPTH=0.    
        THARG=B*(THM-THU)   
        IF(THARG.GT.-20.) EXPTH=EXP(THARG)  
  150   TH=THU+LOG(EXPTH+(1.-EXPTH)*RLU_HIJING(0))/B   
        TH=MAX(THM,MIN(THU,TH)) 
        RATLOG=MIN((B+C*(TH+THM))*(TH-THM),(B+C*(TH+THU))*(TH-THU)) 
        IF(RATLOG.LT.LOG(RLU_HIJING(0))) GOTO 150  
        VINT(21)=1. 
        VINT(22)=0. 
        VINT(23)=MIN(1.,MAX(-1.,(2.*TH-THTER1)/THTER2)) 
    
C...Note: in the following, by In is meant the integral over the    
C...quantity multiplying coefficient cn.    
C...Choose tau according to h1(tau)/tau, where  
C...h1(tau) = c0 + I0/I1*c1*1/tau + I0/I2*c2*1/(tau+tau_R) +    
C...I0/I3*c3*tau/((s*tau-m^2)^2+(m*Gamma)^2) +  
C...I0/I4*c4*1/(tau+tau_R') +   
C...I0/I5*c5*tau/((s*tau-m'^2)^2+(m'*Gamma')^2), and    
C...c0 + c1 + c2 + c3 + c4 + c5 = 1 
      ELSEIF(ISET(ISUB).GE.1.AND.ISET(ISUB).LE.4) THEN  
        CALL PYKLIM_HIJING(1)  
        IF(MINT(51).NE.0) GOTO 100  
        RTAU=RLU_HIJING(0) 
        MTAU=1  
        IF(RTAU.GT.COEF(ISUB,1)) MTAU=2 
        IF(RTAU.GT.COEF(ISUB,1)+COEF(ISUB,2)) MTAU=3    
        IF(RTAU.GT.COEF(ISUB,1)+COEF(ISUB,2)+COEF(ISUB,3)) MTAU=4   
        IF(RTAU.GT.COEF(ISUB,1)+COEF(ISUB,2)+COEF(ISUB,3)+COEF(ISUB,4)) 
     &  MTAU=5  
        IF(RTAU.GT.COEF(ISUB,1)+COEF(ISUB,2)+COEF(ISUB,3)+COEF(ISUB,4)+ 
     &  COEF(ISUB,5)) MTAU=6    
        CALL PYKMAP_HIJING(1,MTAU,RLU_HIJING(0))  
    
C...2 -> 3, 4 processes:    
C...Choose tau' according to h4(tau,tau')/tau', where   
C...h4(tau,tau') = c0 + I0/I1*c1*(1 - tau/tau')^3/tau', and 
C...c0 + c1 = 1.    
        IF(ISET(ISUB).EQ.3.OR.ISET(ISUB).EQ.4) THEN 
          CALL PYKLIM_HIJING(4)    
          IF(MINT(51).NE.0) GOTO 100    
          RTAUP=RLU_HIJING(0)  
          MTAUP=1   
          IF(RTAUP.GT.COEF(ISUB,15)) MTAUP=2    
          CALL PYKMAP_HIJING(4,MTAUP,RLU_HIJING(0))   
        ENDIF   
    
C...Choose y* according to h2(y*), where    
C...h2(y*) = I0/I1*c1*(y*-y*min) + I0/I2*c2*(y*max-y*) +    
C...I0/I3*c3*1/cosh(y*), I0 = y*max-y*min, and c1 + c2 + c3 = 1.    
        CALL PYKLIM_HIJING(2)  
        IF(MINT(51).NE.0) GOTO 100  
        RYST=RLU_HIJING(0) 
        MYST=1  
        IF(RYST.GT.COEF(ISUB,7)) MYST=2 
        IF(RYST.GT.COEF(ISUB,7)+COEF(ISUB,8)) MYST=3    
        CALL PYKMAP_HIJING(2,MYST,RLU_HIJING(0))  
    
C...2 -> 2 processes:   
C...Choose cos(theta-hat) (cth) according to h3(cth), where 
C...h3(cth) = c0 + I0/I1*c1*1/(A - cth) + I0/I2*c2*1/(A + cth) +    
C...I0/I3*c3*1/(A - cth)^2 + I0/I4*c4*1/(A + cth)^2,    
C...A = 1 + 2*(m3*m4/sh)^2 (= 1 for massless products), 
C...and c0 + c1 + c2 + c3 + c4 = 1. 
        CALL PYKLIM_HIJING(3)  
        IF(MINT(51).NE.0) GOTO 100  
        IF(ISET(ISUB).EQ.2.OR.ISET(ISUB).EQ.4) THEN 
          RCTH=RLU_HIJING(0)   
          MCTH=1    
          IF(RCTH.GT.COEF(ISUB,10)) MCTH=2  
          IF(RCTH.GT.COEF(ISUB,10)+COEF(ISUB,11)) MCTH=3    
          IF(RCTH.GT.COEF(ISUB,10)+COEF(ISUB,11)+COEF(ISUB,12)) MCTH=4  
          IF(RCTH.GT.COEF(ISUB,10)+COEF(ISUB,11)+COEF(ISUB,12)+ 
     &    COEF(ISUB,13)) MCTH=5 
          CALL PYKMAP_HIJING(3,MCTH,RLU_HIJING(0))    
        ENDIF   
    
C...Low-pT or multiple interactions (first semihard interaction).   
      ELSEIF(ISET(ISUB).EQ.5) THEN  
        CALL PYMULT_HIJING(3)  
        ISUB=MINT(1)    
      ENDIF 
    
C...Choose azimuthal angle. 
      VINT(24)=PARU(2)*RLU_HIJING(0)   
    
C...Check against user cuts on kinematics at parton level.  
      MINT(51)=0    
      IF(ISUB.LE.90.OR.ISUB.GT.100) CALL PYKLIM_HIJING(0)  
      IF(MINT(51).NE.0) GOTO 100    
      IF(MINT(82).EQ.1.AND.MSTP(141).GE.1) THEN 
        MCUT=0  
        IF(MSUB(91)+MSUB(92)+MSUB(93)+MSUB(94)+MSUB(95).EQ.0)   
     &  CALL PYKCUT_HIJING(MCUT)   
        IF(MCUT.NE.0) GOTO 100  
      ENDIF 
    
C...Calculate differential cross-section for different subprocesses.    
      CALL PYSIGH_HIJING(NCHN,SIGS)    
    
C...Calculations for Monte Carlo estimate of all cross-sections.    
      IF(MINT(82).EQ.1.AND.ISUB.LE.90.OR.ISUB.GE.96) THEN   
        XSEC(ISUB,2)=XSEC(ISUB,2)+SIGS  
      ELSEIF(MINT(82).EQ.1) THEN    
        XSEC(ISUB,2)=XSEC(ISUB,2)+XSEC(ISUB,1)  
      ENDIF 
    
C...Multiple interactions: store results of cross-section calculation.  
      IF(MINT(43).EQ.4.AND.MSTP(82).GE.3) THEN  
        VINT(153)=SIGS  
        CALL PYMULT_HIJING(4)  
      ENDIF 
    
C...Weighting using estimate of maximum of differential cross-section.  
      VIOL=SIGS/XSEC(ISUB,1)    
      IF(VIOL.LT.RLU_HIJING(0)) GOTO 100   
    
C...Check for possible violation of estimated maximum of differential   
C...cross-section used in weighting.    
      IF(MSTP(123).LE.0) THEN   
        IF(VIOL.GT.1.) THEN 
          WRITE(MSTU(11),1000) VIOL,NGEN(0,3)+1 
          WRITE(MSTU(11),1100) ISUB,VINT(21),VINT(22),VINT(23),VINT(26) 
          STOP  
        ENDIF   
      ELSEIF(MSTP(123).EQ.1) THEN   
        IF(VIOL.GT.VINT(108)) THEN  
          VINT(108)=VIOL    
C          IF(VIOL.GT.1.) THEN   
C            WRITE(MSTU(11),1200) VIOL,NGEN(0,3)+1   
C            WRITE(MSTU(11),1100) ISUB,VINT(21),VINT(22),VINT(23),   
C     &      VINT(26)    
C          ENDIF 
        ENDIF   
      ELSEIF(VIOL.GT.VINT(108)) THEN    
        VINT(108)=VIOL  
        IF(VIOL.GT.1.) THEN 
          XDIF=XSEC(ISUB,1)*(VIOL-1.)   
          XSEC(ISUB,1)=XSEC(ISUB,1)+XDIF    
          IF(MSUB(ISUB).EQ.1.AND.(ISUB.LE.90.OR.ISUB.GT.96))    
     &    XSEC(0,1)=XSEC(0,1)+XDIF  
C          WRITE(MSTU(11),1200) VIOL,NGEN(0,3)+1 
C          WRITE(MSTU(11),1100) ISUB,VINT(21),VINT(22),VINT(23),VINT(26) 
C          IF(ISUB.LE.9) THEN    
C            WRITE(MSTU(11),1300) ISUB,XSEC(ISUB,1)  
C          ELSEIF(ISUB.LE.99) THEN   
C            WRITE(MSTU(11),1400) ISUB,XSEC(ISUB,1)  
C          ELSE  
C            WRITE(MSTU(11),1500) ISUB,XSEC(ISUB,1)  
C          ENDIF 
          VINT(108)=1.  
        ENDIF   
      ENDIF 
    
C...Multiple interactions: choose impact parameter. 
      VINT(148)=1.  
      IF(MINT(43).EQ.4.AND.(ISUB.LE.90.OR.ISUB.GE.96).AND.MSTP(82).GE.3)    
     &THEN  
        CALL PYMULT_HIJING(5)  
        IF(VINT(150).LT.RLU_HIJING(0)) GOTO 100    
      ENDIF 
      IF(MINT(82).EQ.1.AND.MSUB(95).EQ.1) THEN  
        IF(ISUB.LE.90.OR.ISUB.GE.95) NGEN(95,1)=NGEN(95,1)+1    
        IF(ISUB.LE.90.OR.ISUB.GE.96) NGEN(96,2)=NGEN(96,2)+1    
      ENDIF 
      IF(ISUB.LE.90.OR.ISUB.GE.96) MINT(31)=MINT(31)+1  
    
C...Choose flavour of reacting partons (and subprocess).    
      RSIGS=SIGS*RLU_HIJING(0) 
      QT2=VINT(48)  
      RQQBAR=PARP(87)*(1.-(QT2/(QT2+(PARP(88)*PARP(82))**2))**2)    
      IF(ISUB.NE.95.AND.(ISUB.NE.96.OR.MSTP(82).LE.1.OR.    
     &RLU_HIJING(0).GT.RQQBAR)) THEN   
        DO 190 ICHN=1,NCHN  
        KFL1=ISIG(ICHN,1)   
        KFL2=ISIG(ICHN,2)   
        MINT(2)=ISIG(ICHN,3)    
        RSIGS=RSIGS-SIGH(ICHN)  
        IF(RSIGS.LE.0.) GOTO 210    
  190   CONTINUE    
    
C...Multiple interactions: choose qqbar preferentially at small pT. 
      ELSEIF(ISUB.EQ.96) THEN   
        CALL PYSPLI_HIJING(MINT(11),21,KFL1,KFLDUM)    
        CALL PYSPLI_HIJING(MINT(12),21,KFL2,KFLDUM)    
        MINT(1)=11  
        MINT(2)=1   
        IF(KFL1.EQ.KFL2.AND.RLU_HIJING(0).LT.0.5) MINT(2)=2    
    
C...Low-pT: choose string drawing configuration.    
      ELSE  
        KFL1=21 
        KFL2=21 
        RSIGS=6.*RLU_HIJING(0) 
        MINT(2)=1   
        IF(RSIGS.GT.1.) MINT(2)=2   
        IF(RSIGS.GT.2.) MINT(2)=3   
      ENDIF 
    
C...Reassign QCD process. Partons before initial state radiation.   
  210 IF(MINT(2).GT.10) THEN    
        MINT(1)=MINT(2)/10  
        MINT(2)=MOD(MINT(2),10) 
      ENDIF 
      MINT(15)=KFL1 
      MINT(16)=KFL2 
      MINT(13)=MINT(15) 
      MINT(14)=MINT(16) 
      VINT(141)=VINT(41)    
      VINT(142)=VINT(42)    
    
C...Format statements for differential cross-section maximum violations.    
 1000 FORMAT(1X,'Error: maximum violated by',1P,E11.3,1X,   
     &'in event',1X,I7,'.'/1X,'Execution stopped!') 
 1100 FORMAT(1X,'ISUB = ',I3,'; Point of violation:'/1X,'tau =',1P, 
     &E11.3,', y* =',E11.3,', cthe = ',0P,F11.7,', tau'' =',1P,E11.3)   
 1200 FORMAT(1X,'Warning: maximum violated by',1P,E11.3,1X, 
     &'in event',1X,I7) 
 1300 FORMAT(1X,'XSEC(',I1,',1) increased to',1P,E11.3) 
 1400 FORMAT(1X,'XSEC(',I2,',1) increased to',1P,E11.3) 
 1500 FORMAT(1X,'XSEC(',I3,',1) increased to',1P,E11.3) 
    
      RETURN    
      END   
