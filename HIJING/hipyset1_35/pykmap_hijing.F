* $Id$
    
C*********************************************************************  
    
      SUBROUTINE PYKMAP_HIJING(IVAR,MVAR,VVAR) 
    
C...Maps a uniform distribution into a distribution of a kinematical    
C...variable according to one of the possibilities allowed. It is   
C...assumed that kinematical limits have been set by a PYKLIM call. 
#include "ludat2_hijing.inc"
#include "pyint1_hijing.inc"
#include "pyint2_hijing.inc"
    
C...Convert VVAR to tau variable.   
      ISUB=MINT(1)  
      IF(IVAR.EQ.1) THEN    
        TAUMIN=VINT(11) 
        TAUMAX=VINT(31) 
        IF(MVAR.EQ.3.OR.MVAR.EQ.4) THEN 
          TAURE=VINT(73)    
          GAMRE=VINT(74)    
        ELSEIF(MVAR.EQ.5.OR.MVAR.EQ.6) THEN 
          TAURE=VINT(75)    
          GAMRE=VINT(76)    
        ENDIF   
        IF(MINT(43).EQ.1.AND.(ISET(ISUB).EQ.1.OR.ISET(ISUB).EQ.2)) THEN 
          TAU=1.    
        ELSEIF(MVAR.EQ.1) THEN  
          TAU=TAUMIN*(TAUMAX/TAUMIN)**VVAR  
        ELSEIF(MVAR.EQ.2) THEN  
          TAU=TAUMAX*TAUMIN/(TAUMIN+(TAUMAX-TAUMIN)*VVAR)   
        ELSEIF(MVAR.EQ.3.OR.MVAR.EQ.5) THEN 
          RATGEN=(TAURE+TAUMAX)/(TAURE+TAUMIN)*TAUMIN/TAUMAX    
          TAU=TAURE*TAUMIN/((TAURE+TAUMIN)*RATGEN**VVAR-TAUMIN) 
        ELSE    
          AUPP=ATAN((TAUMAX-TAURE)/GAMRE)   
          ALOW=ATAN((TAUMIN-TAURE)/GAMRE)   
          TAU=TAURE+GAMRE*TAN(ALOW+(AUPP-ALOW)*VVAR)    
        ENDIF   
        VINT(21)=MIN(TAUMAX,MAX(TAUMIN,TAU))    
    
C...Convert VVAR to y* variable.    
      ELSEIF(IVAR.EQ.2) THEN    
        YSTMIN=VINT(12) 
        YSTMAX=VINT(32) 
        IF(MINT(43).EQ.1) THEN  
          YST=0.    
        ELSEIF(MINT(43).EQ.2) THEN  
          IF(ISET(ISUB).LE.2) YST=-0.5*LOG(VINT(21))    
          IF(ISET(ISUB).GE.3) YST=-0.5*LOG(VINT(26))    
        ELSEIF(MINT(43).EQ.3) THEN  
          IF(ISET(ISUB).LE.2) YST=0.5*LOG(VINT(21)) 
          IF(ISET(ISUB).GE.3) YST=0.5*LOG(VINT(26)) 
        ELSEIF(MVAR.EQ.1) THEN  
          YST=YSTMIN+(YSTMAX-YSTMIN)*SQRT(VVAR) 
        ELSEIF(MVAR.EQ.2) THEN  
          YST=YSTMAX-(YSTMAX-YSTMIN)*SQRT(1.-VVAR)  
        ELSE    
          AUPP=ATAN(EXP(YSTMAX))    
          ALOW=ATAN(EXP(YSTMIN))    
          YST=LOG(TAN(ALOW+(AUPP-ALOW)*VVAR))   
        ENDIF   
        VINT(22)=MIN(YSTMAX,MAX(YSTMIN,YST))    
    
C...Convert VVAR to cos(theta-hat) variable.    
      ELSEIF(IVAR.EQ.3) THEN    
        RM34=2.*VINT(63)*VINT(64)/(VINT(21)*VINT(2))**2 
        RSQM=1.+RM34    
        IF(2.*VINT(71)**2/(VINT(21)*VINT(2)).LT.0.0001) RM34=MAX(RM34,  
     &  2.*VINT(71)**2/(VINT(21)*VINT(2)))  
        CTNMIN=VINT(13) 
        CTNMAX=VINT(33) 
        CTPMIN=VINT(14) 
        CTPMAX=VINT(34) 
        IF(MVAR.EQ.1) THEN  
          ANEG=CTNMAX-CTNMIN    
          APOS=CTPMAX-CTPMIN    
          IF(ANEG.GT.0..AND.VVAR*(ANEG+APOS).LE.ANEG) THEN  
            VCTN=VVAR*(ANEG+APOS)/ANEG  
            CTH=CTNMIN+(CTNMAX-CTNMIN)*VCTN 
          ELSE  
            VCTP=(VVAR*(ANEG+APOS)-ANEG)/APOS   
            CTH=CTPMIN+(CTPMAX-CTPMIN)*VCTP 
          ENDIF 
        ELSEIF(MVAR.EQ.2) THEN  
          RMNMIN=MAX(RM34,RSQM-CTNMIN)  
          RMNMAX=MAX(RM34,RSQM-CTNMAX)  
          RMPMIN=MAX(RM34,RSQM-CTPMIN)  
          RMPMAX=MAX(RM34,RSQM-CTPMAX)  
          ANEG=LOG(RMNMIN/RMNMAX)   
          APOS=LOG(RMPMIN/RMPMAX)   
          IF(ANEG.GT.0..AND.VVAR*(ANEG+APOS).LE.ANEG) THEN  
            VCTN=VVAR*(ANEG+APOS)/ANEG  
            CTH=RSQM-RMNMIN*(RMNMAX/RMNMIN)**VCTN   
          ELSE  
            VCTP=(VVAR*(ANEG+APOS)-ANEG)/APOS   
            CTH=RSQM-RMPMIN*(RMPMAX/RMPMIN)**VCTP   
          ENDIF 
        ELSEIF(MVAR.EQ.3) THEN  
          RMNMIN=MAX(RM34,RSQM+CTNMIN)  
          RMNMAX=MAX(RM34,RSQM+CTNMAX)  
          RMPMIN=MAX(RM34,RSQM+CTPMIN)  
          RMPMAX=MAX(RM34,RSQM+CTPMAX)  
          ANEG=LOG(RMNMAX/RMNMIN)   
          APOS=LOG(RMPMAX/RMPMIN)   
          IF(ANEG.GT.0..AND.VVAR*(ANEG+APOS).LE.ANEG) THEN  
            VCTN=VVAR*(ANEG+APOS)/ANEG  
            CTH=RMNMIN*(RMNMAX/RMNMIN)**VCTN-RSQM   
          ELSE  
            VCTP=(VVAR*(ANEG+APOS)-ANEG)/APOS   
            CTH=RMPMIN*(RMPMAX/RMPMIN)**VCTP-RSQM   
          ENDIF 
        ELSEIF(MVAR.EQ.4) THEN  
          RMNMIN=MAX(RM34,RSQM-CTNMIN)  
          RMNMAX=MAX(RM34,RSQM-CTNMAX)  
          RMPMIN=MAX(RM34,RSQM-CTPMIN)  
          RMPMAX=MAX(RM34,RSQM-CTPMAX)  
          ANEG=1./RMNMAX-1./RMNMIN  
          APOS=1./RMPMAX-1./RMPMIN  
          IF(ANEG.GT.0..AND.VVAR*(ANEG+APOS).LE.ANEG) THEN  
            VCTN=VVAR*(ANEG+APOS)/ANEG  
            CTH=RSQM-1./(1./RMNMIN+ANEG*VCTN)   
          ELSE  
            VCTP=(VVAR*(ANEG+APOS)-ANEG)/APOS   
            CTH=RSQM-1./(1./RMPMIN+APOS*VCTP)   
          ENDIF 
        ELSEIF(MVAR.EQ.5) THEN  
          RMNMIN=MAX(RM34,RSQM+CTNMIN)  
          RMNMAX=MAX(RM34,RSQM+CTNMAX)  
          RMPMIN=MAX(RM34,RSQM+CTPMIN)  
          RMPMAX=MAX(RM34,RSQM+CTPMAX)  
          ANEG=1./RMNMIN-1./RMNMAX  
          APOS=1./RMPMIN-1./RMPMAX  
          IF(ANEG.GT.0..AND.VVAR*(ANEG+APOS).LE.ANEG) THEN  
            VCTN=VVAR*(ANEG+APOS)/ANEG  
            CTH=1./(1./RMNMIN-ANEG*VCTN)-RSQM   
          ELSE  
            VCTP=(VVAR*(ANEG+APOS)-ANEG)/APOS   
            CTH=1./(1./RMPMIN-APOS*VCTP)-RSQM   
          ENDIF 
        ENDIF   
        IF(CTH.LT.0.) CTH=MIN(CTNMAX,MAX(CTNMIN,CTH))   
        IF(CTH.GT.0.) CTH=MIN(CTPMAX,MAX(CTPMIN,CTH))   
        VINT(23)=CTH    
    
C...Convert VVAR to tau' variable.  
      ELSEIF(IVAR.EQ.4) THEN    
        TAU=VINT(11)    
        TAUPMN=VINT(16) 
        TAUPMX=VINT(36) 
        IF(MINT(43).EQ.1) THEN  
          TAUP=1.   
        ELSEIF(MVAR.EQ.1) THEN  
          TAUP=TAUPMN*(TAUPMX/TAUPMN)**VVAR 
        ELSE    
          AUPP=(1.-TAU/TAUPMX)**4   
          ALOW=(1.-TAU/TAUPMN)**4   
          TAUP=TAU/(1.-(ALOW+(AUPP-ALOW)*VVAR)**0.25)   
        ENDIF   
        VINT(26)=MIN(TAUPMX,MAX(TAUPMN,TAUP))   
      ENDIF 
    
      RETURN    
      END   
