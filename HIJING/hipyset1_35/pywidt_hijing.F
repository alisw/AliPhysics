* $Id$
    
C*********************************************************************  
    
      SUBROUTINE PYWIDT_HIJING(KFLR,RMAS,WDTP,WDTE)    
    
C...Calculates full and partial widths of resonances.   
#include "ludat1_hijing.inc"
#include "ludat2_hijing.inc"
#include "ludat3_hijing.inc"
#include "pypars_hijing.inc"
#include "pyint1_hijing.inc"
#include "pyint4_hijing.inc"
      DIMENSION WDTP(0:40),WDTE(0:40,0:5)   
    
C...Some common constants.  
      KFLA=IABS(KFLR)   
      SQM=RMAS**2   
      AS=ULALPS_HIJING(SQM)    
      AEM=PARU(101) 
      XW=PARU(102)  
      RADC=1.+AS/PARU(1)    
    
C...Reset width information.    
      DO 100 I=0,40 
      WDTP(I)=0.    
      DO 100 J=0,5  
  100 WDTE(I,J)=0.  
    
      IF(KFLA.EQ.21) THEN   
C...QCD:    
        DO 110 I=1,MDCY(21,3)   
        IDC=I+MDCY(21,2)-1  
        RM1=(PMAS(IABS(KFDP(IDC,1)),1)/RMAS)**2 
        RM2=(PMAS(IABS(KFDP(IDC,2)),1)/RMAS)**2 
        IF(SQRT(RM1)+SQRT(RM2).GT.1..OR.MDME(IDC,1).LT.0) GOTO 110  
        IF(I.LE.8) THEN 
C...QCD -> q + qb   
          WDTP(I)=(1.+2.*RM1)*SQRT(MAX(0.,1.-4.*RM1))   
          WID2=1.   
        ENDIF   
        WDTP(0)=WDTP(0)+WDTP(I) 
        IF(MDME(IDC,1).GT.0) THEN   
          WDTE(I,MDME(IDC,1))=WDTP(I)*WID2  
          WDTE(0,MDME(IDC,1))=WDTE(0,MDME(IDC,1))+WDTE(I,MDME(IDC,1))   
          WDTE(I,0)=WDTE(I,MDME(IDC,1)) 
          WDTE(0,0)=WDTE(0,0)+WDTE(I,0) 
        ENDIF   
  110   CONTINUE    
    
      ELSEIF(KFLA.EQ.23) THEN   
C...Z0: 
        IF(MINT(61).EQ.1) THEN  
          EI=KCHG(IABS(MINT(15)),1)/3.  
          AI=SIGN(1.,EI)    
          VI=AI-4.*EI*XW    
          SQMZ=PMAS(23,1)**2    
          GZMZ=PMAS(23,2)*PMAS(23,1)    
          GGI=EI**2 
          GZI=EI*VI/(8.*XW*(1.-XW))*SQM*(SQM-SQMZ)/ 
     &    ((SQM-SQMZ)**2+GZMZ**2)   
          ZZI=(VI**2+AI**2)/(16.*XW*(1.-XW))**2*SQM**2/ 
     &    ((SQM-SQMZ)**2+GZMZ**2)   
          IF(MSTP(43).EQ.1) THEN    
C...Only gamma* production included 
            GZI=0.  
            ZZI=0.  
          ELSEIF(MSTP(43).EQ.2) THEN    
C...Only Z0 production included 
            GGI=0.  
            GZI=0.  
          ENDIF 
        ELSEIF(MINT(61).EQ.2) THEN  
          VINT(111)=0.  
          VINT(112)=0.  
          VINT(114)=0.  
        ENDIF   
        DO 120 I=1,MDCY(23,3)   
        IDC=I+MDCY(23,2)-1  
        RM1=(PMAS(IABS(KFDP(IDC,1)),1)/RMAS)**2 
        RM2=(PMAS(IABS(KFDP(IDC,2)),1)/RMAS)**2 
        IF(SQRT(RM1)+SQRT(RM2).GT.1..OR.MDME(IDC,1).LT.0) GOTO 120  
        IF(I.LE.8) THEN 
C...Z0 -> q + qb    
          EF=KCHG(I,1)/3.   
          AF=SIGN(1.,EF+0.1)    
          VF=AF-4.*EF*XW    
          IF(MINT(61).EQ.0) THEN    
            WDTP(I)=3.*(VF**2*(1.+2.*RM1)+AF**2*(1.-4.*RM1))*   
     &      SQRT(MAX(0.,1.-4.*RM1))*RADC    
          ELSEIF(MINT(61).EQ.1) THEN    
            WDTP(I)=3.*((GGI*EF**2+GZI*EF*VF+ZZI*VF**2)*    
     &      (1.+2.*RM1)+ZZI*AF**2*(1.-4.*RM1))* 
     &      SQRT(MAX(0.,1.-4.*RM1))*RADC    
          ELSEIF(MINT(61).EQ.2) THEN    
            GGF=3.*EF**2*(1.+2.*RM1)*SQRT(MAX(0.,1.-4.*RM1))*RADC   
            GZF=3.*EF*VF*(1.+2.*RM1)*SQRT(MAX(0.,1.-4.*RM1))*RADC   
            ZZF=3.*(VF**2*(1.+2.*RM1)+AF**2*(1.-4.*RM1))*   
     &      SQRT(MAX(0.,1.-4.*RM1))*RADC    
          ENDIF 
          WID2=1.   
        ELSEIF(I.LE.16) THEN    
C...Z0 -> l+ + l-, nu + nub 
          EF=KCHG(I+2,1)/3. 
          AF=SIGN(1.,EF+0.1)    
          VF=AF-4.*EF*XW    
          WDTP(I)=(VF**2*(1.+2.*RM1)+AF**2*(1.-4.*RM1))*    
     &    SQRT(MAX(0.,1.-4.*RM1))   
          IF(MINT(61).EQ.0) THEN    
            WDTP(I)=(VF**2*(1.+2.*RM1)+AF**2*(1.-4.*RM1))*  
     &      SQRT(MAX(0.,1.-4.*RM1)) 
          ELSEIF(MINT(61).EQ.1) THEN    
            WDTP(I)=((GGI*EF**2+GZI*EF*VF+ZZI*VF**2)*   
     &      (1.+2.*RM1)+ZZI*AF**2*(1.-4.*RM1))* 
     &      SQRT(MAX(0.,1.-4.*RM1)) 
          ELSEIF(MINT(61).EQ.2) THEN    
            GGF=EF**2*(1.+2.*RM1)*SQRT(MAX(0.,1.-4.*RM1))   
            GZF=EF*VF*(1.+2.*RM1)*SQRT(MAX(0.,1.-4.*RM1))   
            ZZF=(VF**2*(1.+2.*RM1)+AF**2*(1.-4.*RM1))*  
     &      SQRT(MAX(0.,1.-4.*RM1)) 
          ENDIF 
          WID2=1.   
        ELSE    
C...Z0 -> H+ + H-   
          CF=2.*(1.-2.*XW)  
          IF(MINT(61).EQ.0) THEN    
            WDTP(I)=0.25*CF**2*(1.-4.*RM1)*SQRT(MAX(0.,1.-4.*RM1))  
          ELSEIF(MINT(61).EQ.1) THEN    
            WDTP(I)=0.25*(GGI+GZI*CF+ZZI*CF**2)*(1.-4.*RM1)*    
     &      SQRT(MAX(0.,1.-4.*RM1)) 
          ELSEIF(MINT(61).EQ.2) THEN    
            GGF=0.25*(1.-4.*RM1)*SQRT(MAX(0.,1.-4.*RM1))    
            GZF=0.25*CF*(1.-4.*RM1)*SQRT(MAX(0.,1.-4.*RM1)) 
            ZZF=0.25*CF**2*(1.-4.*RM1)*SQRT(MAX(0.,1.-4.*RM1))  
          ENDIF 
          WID2=WIDS(37,1)   
        ENDIF   
        WDTP(0)=WDTP(0)+WDTP(I) 
        IF(MDME(IDC,1).GT.0) THEN   
          WDTE(I,MDME(IDC,1))=WDTP(I)*WID2  
          WDTE(0,MDME(IDC,1))=WDTE(0,MDME(IDC,1))+WDTE(I,MDME(IDC,1))   
          WDTE(I,0)=WDTE(I,MDME(IDC,1)) 
          WDTE(0,0)=WDTE(0,0)+WDTE(I,0) 
          VINT(111)=VINT(111)+GGF*WID2  
          VINT(112)=VINT(112)+GZF*WID2  
          VINT(114)=VINT(114)+ZZF*WID2  
        ENDIF   
  120   CONTINUE    
        IF(MSTP(43).EQ.1) THEN  
C...Only gamma* production included 
          VINT(112)=0.  
          VINT(114)=0.  
        ELSEIF(MSTP(43).EQ.2) THEN  
C...Only Z0 production included 
          VINT(111)=0.  
          VINT(112)=0.  
        ENDIF   
    
      ELSEIF(KFLA.EQ.24) THEN   
C...W+/-:   
        DO 130 I=1,MDCY(24,3)   
        IDC=I+MDCY(24,2)-1  
        RM1=(PMAS(IABS(KFDP(IDC,1)),1)/RMAS)**2 
        RM2=(PMAS(IABS(KFDP(IDC,2)),1)/RMAS)**2 
        IF(SQRT(RM1)+SQRT(RM2).GT.1..OR.MDME(IDC,1).LT.0) GOTO 130  
        IF(I.LE.16) THEN    
C...W+/- -> q + qb' 
          WDTP(I)=3.*(2.-RM1-RM2-(RM1-RM2)**2)* 
     &    SQRT(MAX(0.,(1.-RM1-RM2)**2-4.*RM1*RM2))* 
     &    VCKM((I-1)/4+1,MOD(I-1,4)+1)*RADC 
          WID2=1.   
        ELSE    
C...W+/- -> l+/- + nu   
          WDTP(I)=(2.-RM1-RM2-(RM1-RM2)**2)*    
     &    SQRT(MAX(0.,(1.-RM1-RM2)**2-4.*RM1*RM2))  
          WID2=1.   
        ENDIF   
        WDTP(0)=WDTP(0)+WDTP(I) 
        IF(MDME(IDC,1).GT.0) THEN   
          WDTE(I,MDME(IDC,1))=WDTP(I)*WID2  
          WDTE(0,MDME(IDC,1))=WDTE(0,MDME(IDC,1))+WDTE(I,MDME(IDC,1))   
          WDTE(I,0)=WDTE(I,MDME(IDC,1)) 
          WDTE(0,0)=WDTE(0,0)+WDTE(I,0) 
        ENDIF   
  130   CONTINUE    
    
      ELSEIF(KFLA.EQ.25) THEN   
C...H0: 
        DO 170 I=1,MDCY(25,3)   
        IDC=I+MDCY(25,2)-1  
        RM1=(PMAS(IABS(KFDP(IDC,1)),1)/RMAS)**2 
        RM2=(PMAS(IABS(KFDP(IDC,2)),1)/RMAS)**2 
        IF(SQRT(RM1)+SQRT(RM2).GT.1..OR.MDME(IDC,1).LT.0) GOTO 170  
        IF(I.LE.8) THEN 
C...H0 -> q + qb    
          WDTP(I)=3.*RM1*(1.-4.*RM1)*SQRT(MAX(0.,1.-4.*RM1))*RADC   
          WID2=1.   
        ELSEIF(I.LE.12) THEN    
C...H0 -> l+ + l-   
          WDTP(I)=RM1*(1.-4.*RM1)*SQRT(MAX(0.,1.-4.*RM1))   
          WID2=1.   
        ELSEIF(I.EQ.13) THEN    
C...H0 -> g + g; quark loop contribution only   
          ETARE=0.  
          ETAIM=0.  
          DO 140 J=1,2*MSTP(1)  
          EPS=(2.*PMAS(J,1)/RMAS)**2    
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
  140     CONTINUE  
          ETA2=ETARE**2+ETAIM**2    
          WDTP(I)=(AS/PARU(1))**2*ETA2  
          WID2=1.   
        ELSEIF(I.EQ.14) THEN    
C...H0 -> gamma + gamma; quark, charged lepton and W loop contributions 
          ETARE=0.  
          ETAIM=0.  
          DO 150 J=1,3*MSTP(1)+1    
          IF(J.LE.2*MSTP(1)) THEN   
            EJ=KCHG(J,1)/3. 
            EPS=(2.*PMAS(J,1)/RMAS)**2  
          ELSEIF(J.LE.3*MSTP(1)) THEN   
            JL=2*(J-2*MSTP(1))-1    
            EJ=KCHG(10+JL,1)/3. 
            EPS=(2.*PMAS(10+JL,1)/RMAS)**2  
          ELSE  
            EPS=(2.*PMAS(24,1)/RMAS)**2 
          ENDIF 
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
          IF(J.LE.2*MSTP(1)) THEN   
            ETARE=ETARE+0.5*3.*EJ**2*EPS*(1.+(EPS-1.)*PHIRE)    
            ETAIM=ETAIM+0.5*3.*EJ**2*EPS*(EPS-1.)*PHIIM 
          ELSEIF(J.LE.3*MSTP(1)) THEN   
            ETARE=ETARE+0.5*EJ**2*EPS*(1.+(EPS-1.)*PHIRE)   
            ETAIM=ETAIM+0.5*EJ**2*EPS*(EPS-1.)*PHIIM    
          ELSE  
            ETARE=ETARE-0.5-0.75*EPS*(1.+(EPS-2.)*PHIRE)    
            ETAIM=ETAIM+0.75*EPS*(EPS-2.)*PHIIM 
          ENDIF 
  150     CONTINUE  
          ETA2=ETARE**2+ETAIM**2    
          WDTP(I)=(AEM/PARU(1))**2*0.5*ETA2 
          WID2=1.   
        ELSEIF(I.EQ.15) THEN    
C...H0 -> gamma + Z0; quark, charged lepton and W loop contributions    
          ETARE=0.  
          ETAIM=0.  
          DO 160 J=1,3*MSTP(1)+1    
          IF(J.LE.2*MSTP(1)) THEN   
            EJ=KCHG(J,1)/3. 
            AJ=SIGN(1.,EJ+0.1)  
            VJ=AJ-4.*EJ*XW  
            EPS=(2.*PMAS(J,1)/RMAS)**2  
            EPSP=(2.*PMAS(J,1)/PMAS(23,1))**2   
          ELSEIF(J.LE.3*MSTP(1)) THEN   
            JL=2*(J-2*MSTP(1))-1    
            EJ=KCHG(10+JL,1)/3. 
            AJ=SIGN(1.,EJ+0.1)  
            VJ=AI-4.*EJ*XW  
            EPS=(2.*PMAS(10+JL,1)/RMAS)**2  
            EPSP=(2.*PMAS(10+JL,1)/PMAS(23,1))**2   
          ELSE  
            EPS=(2.*PMAS(24,1)/RMAS)**2 
            EPSP=(2.*PMAS(24,1)/PMAS(23,1))**2  
          ENDIF 
          IF(EPS.LE.1.) THEN    
            ROOT=SQRT(1.-EPS)   
            IF(EPS.GT.1.E-4) THEN   
              RLN=LOG((1.+ROOT)/(1.-ROOT))  
            ELSE    
              RLN=LOG(4./EPS-2.)    
            ENDIF   
            PHIRE=0.25*(RLN**2-PARU(1)**2)  
            PHIIM=0.5*PARU(1)*RLN   
            PSIRE=-(1.+0.5*ROOT*RLN)    
            PSIIM=0.5*PARU(1)*ROOT  
          ELSE  
            PHIRE=-(ASIN(1./SQRT(EPS)))**2  
            PHIIM=0.    
            PSIRE=-(1.+SQRT(EPS-1.)*ASIN(1./SQRT(EPS))) 
            PSIIM=0.    
          ENDIF 
          IF(EPSP.LE.1.) THEN   
            ROOT=SQRT(1.-EPSP)  
            IF(EPSP.GT.1.E-4) THEN  
              RLN=LOG((1.+ROOT)/(1.-ROOT))  
            ELSE    
              RLN=LOG(4./EPSP-2.)   
            ENDIF   
            PHIREP=0.25*(RLN**2-PARU(1)**2) 
            PHIIMP=0.5*PARU(1)*RLN  
            PSIREP=-(1.+0.5*ROOT*RLN)   
            PSIIMP=0.5*PARU(1)*ROOT 
          ELSE  
            PHIREP=-(ASIN(1./SQRT(EPSP)))**2    
            PHIIMP=0.   
            PSIREP=-(1.+SQRT(EPSP-1.)*ASIN(1./SQRT(EPSP)))  
            PSIIMP=0.   
          ENDIF 
          FXYRE=EPS*EPSP/(8.*(EPS-EPSP))*(1.-EPS*EPSP/(EPS-EPSP)*(PHIRE-    
     &    PHIREP)+2.*EPS/(EPS-EPSP)*(PSIRE-PSIREP)) 
          FXYIM=EPS*EPSP/(8.*(EPS-EPSP))*(-EPS*EPSP/(EPS-EPSP)*(PHIIM-  
     &    PHIIMP)+2.*EPS/(EPS-EPSP)*(PSIIM-PSIIMP)) 
          F1RE=EPS*EPSP/(2.*(EPS-EPSP))*(PHIRE-PHIREP)  
          F1IM=EPS*EPSP/(2.*(EPS-EPSP))*(PHIIM-PHIIMP)  
          IF(J.LE.2*MSTP(1)) THEN   
            ETARE=ETARE-3.*EJ*VJ*(FXYRE-0.25*F1RE)  
            ETAIM=ETAIM-3.*EJ*VJ*(FXYIM-0.25*F1IM)  
          ELSEIF(J.LE.3*MSTP(1)) THEN   
            ETARE=ETARE-EJ*VJ*(FXYRE-0.25*F1RE) 
            ETAIM=ETAIM-EJ*VJ*(FXYIM-0.25*F1IM) 
          ELSE  
            ETARE=ETARE-SQRT(1.-XW)*(((1.+2./EPS)*XW/SQRT(1.-XW)-   
     &      (5.+2./EPS))*FXYRE+(3.-XW/SQRT(1.-XW))*F1RE)    
            ETAIM=ETAIM-SQRT(1.-XW)*(((1.+2./EPS)*XW/SQRT(1.-XW)-   
     &      (5.+2./EPS))*FXYIM+(3.-XW/SQRT(1.-XW))*F1IM)    
          ENDIF 
  160     CONTINUE  
          ETA2=ETARE**2+ETAIM**2    
          WDTP(I)=(AEM/PARU(1))**2*(1.-(PMAS(23,1)/RMAS)**2)**3/XW*ETA2 
          WID2=WIDS(23,2)   
        ELSE    
C...H0 -> Z0 + Z0, W+ + W-  
          WDTP(I)=(1.-4.*RM1+12.*RM1**2)*SQRT(MAX(0.,1.-4.*RM1))/   
     &    (2.*(18-I))   
          WID2=WIDS(7+I,1)  
        ENDIF   
        WDTP(0)=WDTP(0)+WDTP(I) 
        IF(MDME(IDC,1).GT.0) THEN   
          WDTE(I,MDME(IDC,1))=WDTP(I)*WID2  
          WDTE(0,MDME(IDC,1))=WDTE(0,MDME(IDC,1))+WDTE(I,MDME(IDC,1))   
          WDTE(I,0)=WDTE(I,MDME(IDC,1)) 
          WDTE(0,0)=WDTE(0,0)+WDTE(I,0) 
        ENDIF   
  170   CONTINUE    
    
      ELSEIF(KFLA.EQ.32) THEN   
C...Z'0:    
C....  this line added by A.M.
         API = 0.
C....         
        IF(MINT(61).EQ.1) THEN  
          EI=KCHG(IABS(MINT(15)),1)/3.  
          AI=SIGN(1.,EI)    
          VI=AI-4.*EI*XW    
          SQMZ=PMAS(23,1)**2    
          GZMZ=PMAS(23,2)*PMAS(23,1)    
          API=SIGN(1.,EI)   
          VPI=API-4.*EI*XW  
          SQMZP=PMAS(32,1)**2   
          GZPMZP=PMAS(32,2)*PMAS(32,1)  
          GGI=EI**2 
          GZI=EI*VI/(8.*XW*(1.-XW))*SQM*(SQM-SQMZ)/ 
     &    ((SQM-SQMZ)**2+GZMZ**2)   
          GZPI=EI*VPI/(8.*XW*(1.-XW))*SQM*(SQM-SQMZP)/  
     &    ((SQM-SQMZP)**2+GZPMZP**2)    
          ZZI=(VI**2+AI**2)/(16.*XW*(1.-XW))**2*SQM**2/ 
     &    ((SQM-SQMZ)**2+GZMZ**2)   
          ZZPI=2.*(VI*VPI+AI*API)/(16.*XW*(1.-XW))**2*  
     &    SQM**2*((SQM-SQMZ)*(SQM-SQMZP)+GZMZ*GZPMZP)/  
     &    (((SQM-SQMZ)**2+GZMZ**2)*((SQM-SQMZP)**2+GZPMZP**2))  
          ZPZPI=(VPI**2+API**2)/(16.*XW*(1.-XW))**2*SQM**2/ 
     &    ((SQM-SQMZP)**2+GZPMZP**2)    
          IF(MSTP(44).EQ.1) THEN    
C...Only gamma* production included 
            GZI=0.  
            GZPI=0. 
            ZZI=0.  
            ZZPI=0. 
            ZPZPI=0.    
          ELSEIF(MSTP(44).EQ.2) THEN    
C...Only Z0 production included 
            GGI=0.  
            GZI=0.  
            GZPI=0. 
            ZZPI=0. 
            ZPZPI=0.    
          ELSEIF(MSTP(44).EQ.3) THEN    
C...Only Z'0 production included    
            GGI=0.  
            GZI=0.  
            GZPI=0. 
            ZZI=0.  
            ZZPI=0. 
          ELSEIF(MSTP(44).EQ.4) THEN    
C...Only gamma*/Z0 production included  
            GZPI=0. 
            ZZPI=0. 
            ZPZPI=0.    
          ELSEIF(MSTP(44).EQ.5) THEN    
C...Only gamma*/Z'0 production included 
            GZI=0.  
            ZZI=0.  
            ZZPI=0. 
          ELSEIF(MSTP(44).EQ.6) THEN    
C...Only Z0/Z'0 production included 
            GGI=0.  
            GZI=0.  
            GZPI=0. 
          ENDIF 
        ELSEIF(MINT(61).EQ.2) THEN  
          VINT(111)=0.  
          VINT(112)=0.  
          VINT(113)=0.  
          VINT(114)=0.  
          VINT(115)=0.  
          VINT(116)=0.  
        ENDIF   
        DO 180 I=1,MDCY(32,3)   
        IDC=I+MDCY(32,2)-1  
        RM1=(PMAS(IABS(KFDP(IDC,1)),1)/RMAS)**2 
        RM2=(PMAS(IABS(KFDP(IDC,2)),1)/RMAS)**2 
        IF(SQRT(RM1)+SQRT(RM2).GT.1..OR.MDME(IDC,1).LT.0) GOTO 180  
        IF(I.LE.8) THEN 
C...Z'0 -> q + qb   
          EF=KCHG(I,1)/3.   
          AF=SIGN(1.,EF+0.1)    
          VF=AF-4.*EF*XW    
          APF=SIGN(1.,EF+0.1)   
          VPF=APF-4.*EF*XW  
          IF(MINT(61).EQ.0) THEN    
            WDTP(I)=3.*(VPF**2*(1.+2.*RM1)+APF**2*(1.-4.*RM1))* 
     &      SQRT(MAX(0.,1.-4.*RM1))*RADC    
          ELSEIF(MINT(61).EQ.1) THEN    
            WDTP(I)=3.*((GGI*EF**2+GZI*EF*VF+GZPI*EF*VPF+ZZI*VF**2+ 
     &      ZZPI*VF*VPF+ZPZPI*VPF**2)*(1.+2.*RM1)+(ZZI*AF**2+   
     &      ZZPI*AF*APF+ZPZPI*APF**2)*(1.-4.*RM1))* 
     &      SQRT(MAX(0.,1.-4.*RM1))*RADC    
          ELSEIF(MINT(61).EQ.2) THEN    
            GGF=3.*EF**2*(1.+2.*RM1)*SQRT(MAX(0.,1.-4.*RM1))*RADC   
            GZF=3.*EF*VF*(1.+2.*RM1)*SQRT(MAX(0.,1.-4.*RM1))*RADC   
            GZPF=3.*EF*VPF*(1.+2.*RM1)*SQRT(MAX(0.,1.-4.*RM1))*RADC 
            ZZF=3.*(VF**2*(1.+2.*RM1)+AF**2*(1.-4.*RM1))*   
     &      SQRT(MAX(0.,1.-4.*RM1))*RADC    
            ZZPF=3.*(VF*VPF*(1.+2.*RM1)+AF*APF*(1.-4.*RM1))*    
     &      SQRT(MAX(0.,1.-4.*RM1))*RADC    
            ZPZPF=3.*(VPF**2*(1.+2.*RM1)+APF**2*(1.-4.*RM1))*   
     &      SQRT(MAX(0.,1.-4.*RM1))*RADC    
          ENDIF 
          WID2=1.   
        ELSE    
C...Z'0 -> l+ + l-, nu + nub    
          EF=KCHG(I+2,1)/3. 
          AF=SIGN(1.,EF+0.1)    
          VF=AF-4.*EF*XW    
          APF=SIGN(1.,EF+0.1)   
          VPF=API-4.*EF*XW  
          IF(MINT(61).EQ.0) THEN    
            WDTP(I)=(VPF**2*(1.+2.*RM1)+APF**2*(1.-4.*RM1))*    
     &      SQRT(MAX(0.,1.-4.*RM1)) 
          ELSEIF(MINT(61).EQ.1) THEN    
            WDTP(I)=((GGI*EF**2+GZI*EF*VF+GZPI*EF*VPF+ZZI*VF**2+    
     &      ZZPI*VF*VPF+ZPZPI*VPF**2)*(1.+2.*RM1)+(ZZI*AF**2+   
     &      ZZPI*AF*APF+ZPZPI*APF**2)*(1.-4.*RM1))* 
     &      SQRT(MAX(0.,1.-4.*RM1)) 
          ELSEIF(MINT(61).EQ.2) THEN    
            GGF=EF**2*(1.+2.*RM1)*SQRT(MAX(0.,1.-4.*RM1))   
            GZF=EF*VF*(1.+2.*RM1)*SQRT(MAX(0.,1.-4.*RM1))   
            GZPF=EF*VPF*(1.+2.*RM1)*SQRT(MAX(0.,1.-4.*RM1)) 
            ZZF=(VF**2*(1.+2.*RM1)+AF**2*(1.-4.*RM1))*  
     &      SQRT(MAX(0.,1.-4.*RM1)) 
            ZZPF=(VF*VPF*(1.+2.*RM1)+AF*APF*(1.-4.*RM1))*   
     &      SQRT(MAX(0.,1.-4.*RM1)) 
            ZPZPF=(VPF**2*(1.+2.*RM1)+APF**2*(1.-4.*RM1))*  
     &      SQRT(MAX(0.,1.-4.*RM1)) 
          ENDIF 
          WID2=1.   
        ENDIF   
        WDTP(0)=WDTP(0)+WDTP(I) 
        IF(MDME(IDC,1).GT.0) THEN   
          WDTE(I,MDME(IDC,1))=WDTP(I)*WID2  
          WDTE(0,MDME(IDC,1))=WDTE(0,MDME(IDC,1))+WDTE(I,MDME(IDC,1))   
          WDTE(I,0)=WDTE(I,MDME(IDC,1)) 
          WDTE(0,0)=WDTE(0,0)+WDTE(I,0) 
          VINT(111)=VINT(111)+GGF   
          VINT(112)=VINT(112)+GZF   
          VINT(113)=VINT(113)+GZPF  
          VINT(114)=VINT(114)+ZZF   
          VINT(115)=VINT(115)+ZZPF  
          VINT(116)=VINT(116)+ZPZPF 
        ENDIF   
  180   CONTINUE    
        IF(MSTP(44).EQ.1) THEN  
C...Only gamma* production included 
          VINT(112)=0.  
          VINT(113)=0.  
          VINT(114)=0.  
          VINT(115)=0.  
          VINT(116)=0.  
        ELSEIF(MSTP(44).EQ.2) THEN  
C...Only Z0 production included 
          VINT(111)=0.  
          VINT(112)=0.  
          VINT(113)=0.  
          VINT(115)=0.  
          VINT(116)=0.  
        ELSEIF(MSTP(44).EQ.3) THEN  
C...Only Z'0 production included    
          VINT(111)=0.  
          VINT(112)=0.  
          VINT(113)=0.  
          VINT(114)=0.  
          VINT(115)=0.  
        ELSEIF(MSTP(44).EQ.4) THEN  
C...Only gamma*/Z0 production included  
          VINT(113)=0.  
          VINT(115)=0.  
          VINT(116)=0.  
        ELSEIF(MSTP(44).EQ.5) THEN  
C...Only gamma*/Z'0 production included 
          VINT(112)=0.  
          VINT(114)=0.  
          VINT(115)=0.  
        ELSEIF(MSTP(44).EQ.6) THEN  
C...Only Z0/Z'0 production included 
          VINT(111)=0.  
          VINT(112)=0.  
          VINT(113)=0.  
        ENDIF   
    
      ELSEIF(KFLA.EQ.37) THEN   
C...H+/-:   
        DO 190 I=1,MDCY(37,3)   
        IDC=I+MDCY(37,2)-1  
        RM1=(PMAS(IABS(KFDP(IDC,1)),1)/RMAS)**2 
        RM2=(PMAS(IABS(KFDP(IDC,2)),1)/RMAS)**2 
        IF(SQRT(RM1)+SQRT(RM2).GT.1..OR.MDME(IDC,1).LT.0) GOTO 190  
        IF(I.LE.4) THEN 
C...H+/- -> q + qb' 
          WDTP(I)=3.*((RM1*PARU(121)+RM2/PARU(121))*    
     &    (1.-RM1-RM2)-4.*RM1*RM2)* 
     &    SQRT(MAX(0.,(1.-RM1-RM2)**2-4.*RM1*RM2))*RADC 
          WID2=1.   
        ELSE    
C...H+/- -> l+/- + nu   
          WDTP(I)=((RM1*PARU(121)+RM2/PARU(121))*   
     &    (1.-RM1-RM2)-4.*RM1*RM2)* 
     &    SQRT(MAX(0.,(1.-RM1-RM2)**2-4.*RM1*RM2))  
          WID2=1.   
        ENDIF   
        WDTP(0)=WDTP(0)+WDTP(I) 
        IF(MDME(IDC,1).GT.0) THEN   
          WDTE(I,MDME(IDC,1))=WDTP(I)*WID2  
          WDTE(0,MDME(IDC,1))=WDTE(0,MDME(IDC,1))+WDTE(I,MDME(IDC,1))   
          WDTE(I,0)=WDTE(I,MDME(IDC,1)) 
          WDTE(0,0)=WDTE(0,0)+WDTE(I,0) 
        ENDIF   
  190   CONTINUE    
    
      ELSEIF(KFLA.EQ.40) THEN   
C...R:  
        DO 200 I=1,MDCY(40,3)   
        IDC=I+MDCY(40,2)-1  
        RM1=(PMAS(IABS(KFDP(IDC,1)),1)/RMAS)**2 
        RM2=(PMAS(IABS(KFDP(IDC,2)),1)/RMAS)**2 
        IF(SQRT(RM1)+SQRT(RM2).GT.1..OR.MDME(IDC,1).LT.0) GOTO 200  
        IF(I.LE.4) THEN 
C...R -> q + qb'    
          WDTP(I)=3.*RADC   
          WID2=1.   
        ELSE    
C...R -> l+ + l'-   
          WDTP(I)=1.    
          WID2=1.   
        ENDIF   
        WDTP(0)=WDTP(0)+WDTP(I) 
        IF(MDME(IDC,1).GT.0) THEN   
          WDTE(I,MDME(IDC,1))=WDTP(I)*WID2  
          WDTE(0,MDME(IDC,1))=WDTE(0,MDME(IDC,1))+WDTE(I,MDME(IDC,1))   
          WDTE(I,0)=WDTE(I,MDME(IDC,1)) 
          WDTE(0,0)=WDTE(0,0)+WDTE(I,0) 
        ENDIF   
  200   CONTINUE    
    
      ENDIF 
      MINT(61)=0    
    
      RETURN    
      END   
