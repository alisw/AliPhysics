* $Id$
    
C*********************************************************************  
    
      SUBROUTINE PYINRE_HIJING 
    
C...Calculates full and effective widths of guage bosons, stores masses 
C...and widths, rescales coefficients to be used for resonance  
C...production generation.  
#include "ludat1_hijing.inc"
#include "ludat2_hijing.inc"
#include "ludat3_hijing.inc"
#include "pysubs_hijing.inc"
#include "pypars_hijing.inc"
#include "pyint1_hijing.inc"
#include "pyint2_hijing.inc"
#include "pyint4_hijing.inc"
#include "pyint6_hijing.inc"
      DIMENSION WDTP(0:40),WDTE(0:40,0:5)   
    
C...Calculate full and effective widths of gauge bosons.    
      AEM=PARU(101) 
      XW=PARU(102)  
      DO 100 I=21,40    
      DO 100 J=0,40 
      WIDP(I,J)=0.  
  100 WIDE(I,J)=0.  
    
C...W+/-:   
      WMAS=PMAS(24,1)   
      WFAC=AEM/(24.*XW)*WMAS    
      CALL PYWIDT_HIJING(24,WMAS,WDTP,WDTE)    
      WIDS(24,1)=((WDTE(0,1)+WDTE(0,2))*(WDTE(0,1)+WDTE(0,3))+  
     &(WDTE(0,1)+WDTE(0,2)+WDTE(0,1)+WDTE(0,3))*(WDTE(0,4)+WDTE(0,5))+  
     &2.*WDTE(0,4)*WDTE(0,5))/WDTP(0)**2    
      WIDS(24,2)=(WDTE(0,1)+WDTE(0,2)+WDTE(0,4))/WDTP(0)    
      WIDS(24,3)=(WDTE(0,1)+WDTE(0,3)+WDTE(0,4))/WDTP(0)    
      DO 110 I=0,40 
      WIDP(24,I)=WFAC*WDTP(I)   
  110 WIDE(24,I)=WFAC*WDTE(I,0) 
    
C...H+/-:   
      HCMAS=PMAS(37,1)  
      HCFAC=AEM/(8.*XW)*(HCMAS/WMAS)**2*HCMAS   
      CALL PYWIDT_HIJING(37,HCMAS,WDTP,WDTE)   
      WIDS(37,1)=((WDTE(0,1)+WDTE(0,2))*(WDTE(0,1)+WDTE(0,3))+  
     &(WDTE(0,1)+WDTE(0,2)+WDTE(0,1)+WDTE(0,3))*(WDTE(0,4)+WDTE(0,5))+  
     &2.*WDTE(0,4)*WDTE(0,5))/WDTP(0)**2    
      WIDS(37,2)=(WDTE(0,1)+WDTE(0,2)+WDTE(0,4))/WDTP(0)    
      WIDS(37,3)=(WDTE(0,1)+WDTE(0,3)+WDTE(0,4))/WDTP(0)    
      DO 120 I=0,40 
      WIDP(37,I)=HCFAC*WDTP(I)  
  120 WIDE(37,I)=HCFAC*WDTE(I,0)    
    
C...Z0: 
      ZMAS=PMAS(23,1)   
      ZFAC=AEM/(48.*XW*(1.-XW))*ZMAS    
      CALL PYWIDT_HIJING(23,ZMAS,WDTP,WDTE)    
      WIDS(23,1)=((WDTE(0,1)+WDTE(0,2))**2+ 
     &2.*(WDTE(0,1)+WDTE(0,2))*(WDTE(0,4)+WDTE(0,5))+   
     &2.*WDTE(0,4)*WDTE(0,5))/WDTP(0)**2    
      WIDS(23,2)=(WDTE(0,1)+WDTE(0,2)+WDTE(0,4))/WDTP(0)    
      WIDS(23,3)=0. 
      DO 130 I=0,40 
      WIDP(23,I)=ZFAC*WDTP(I)   
  130 WIDE(23,I)=ZFAC*WDTE(I,0) 
    
C...H0: 
      HMAS=PMAS(25,1)   
      HFAC=AEM/(8.*XW)*(HMAS/WMAS)**2*HMAS  
      CALL PYWIDT_HIJING(25,HMAS,WDTP,WDTE)    
      WIDS(25,1)=((WDTE(0,1)+WDTE(0,2))**2+ 
     &2.*(WDTE(0,1)+WDTE(0,2))*(WDTE(0,4)+WDTE(0,5))+   
     &2.*WDTE(0,4)*WDTE(0,5))/WDTP(0)**2    
      WIDS(25,2)=(WDTE(0,1)+WDTE(0,2)+WDTE(0,4))/WDTP(0)    
      WIDS(25,3)=0. 
      DO 140 I=0,40 
      WIDP(25,I)=HFAC*WDTP(I)   
  140 WIDE(25,I)=HFAC*WDTE(I,0) 
    
C...Z'0:    
      ZPMAS=PMAS(32,1)  
      ZPFAC=AEM/(48.*XW*(1.-XW))*ZPMAS  
      CALL PYWIDT_HIJING(32,ZPMAS,WDTP,WDTE)   
      WIDS(32,1)=((WDTE(0,1)+WDTE(0,2)+WDTE(0,3))**2+   
     &2.*(WDTE(0,1)+WDTE(0,2))*(WDTE(0,4)+WDTE(0,5))+   
     &2.*WDTE(0,4)*WDTE(0,5))/WDTP(0)**2    
      WIDS(32,2)=(WDTE(0,1)+WDTE(0,2)+WDTE(0,4))/WDTP(0)    
      WIDS(32,3)=0. 
      DO 150 I=0,40 
      WIDP(32,I)=ZPFAC*WDTP(I)  
  150 WIDE(32,I)=ZPFAC*WDTE(I,0)    
    
C...R:  
      RMAS=PMAS(40,1)   
      RFAC=0.08*RMAS/((MSTP(1)-1)*(1.+6.*(1.+ULALPS_HIJING(RMAS**2)
     $     /PARU(1)))) 
      CALL PYWIDT_HIJING(40,RMAS,WDTP,WDTE)    
      WIDS(40,1)=((WDTE(0,1)+WDTE(0,2))*(WDTE(0,1)+WDTE(0,3))+  
     &(WDTE(0,1)+WDTE(0,2)+WDTE(0,1)+WDTE(0,3))*(WDTE(0,4)+WDTE(0,5))+  
     &2.*WDTE(0,4)*WDTE(0,5))/WDTP(0)**2    
      WIDS(40,2)=(WDTE(0,1)+WDTE(0,2)+WDTE(0,4))/WDTP(0)    
      WIDS(40,3)=(WDTE(0,1)+WDTE(0,3)+WDTE(0,4))/WDTP(0)    
      DO 160 I=0,40 
      WIDP(40,I)=WFAC*WDTP(I)   
  160 WIDE(40,I)=WFAC*WDTE(I,0) 
    
C...Q:  
      KFLQM=1   
      DO 170 I=1,MIN(8,MDCY(21,3))  
      IDC=I+MDCY(21,2)-1    
      IF(MDME(IDC,1).LE.0) GOTO 170 
      KFLQM=I   
  170 CONTINUE  
      MINT(46)=KFLQM    
      KFPR(81,1)=KFLQM  
      KFPR(81,2)=KFLQM  
      KFPR(82,1)=KFLQM  
      KFPR(82,2)=KFLQM  
    
C...Set resonance widths and branching ratios in JETSET.    
      DO 180 I=1,6  
      IF(I.LE.3) KC=I+22    
      IF(I.EQ.4) KC=32  
      IF(I.EQ.5) KC=37  
      IF(I.EQ.6) KC=40  
      PMAS(KC,2)=WIDP(KC,0) 
      PMAS(KC,3)=MIN(0.9*PMAS(KC,1),10.*PMAS(KC,2)) 
      DO 180 J=1,MDCY(KC,3) 
      IDC=J+MDCY(KC,2)-1    
      BRAT(IDC)=WIDE(KC,J)/WIDE(KC,0)   
  180 CONTINUE  
    
C...Special cases in treatment of gamma*/Z0: redefine process name. 
      IF(MSTP(43).EQ.1) THEN    
        PROC(1)='f + fb -> gamma*'  
      ELSEIF(MSTP(43).EQ.2) THEN    
        PROC(1)='f + fb -> Z0'  
      ELSEIF(MSTP(43).EQ.3) THEN    
        PROC(1)='f + fb -> gamma*/Z0'   
      ENDIF 
    
C...Special cases in treatment of gamma*/Z0/Z'0: redefine process name. 
      IF(MSTP(44).EQ.1) THEN    
        PROC(141)='f + fb -> gamma*'    
      ELSEIF(MSTP(44).EQ.2) THEN    
        PROC(141)='f + fb -> Z0'    
      ELSEIF(MSTP(44).EQ.3) THEN    
        PROC(141)='f + fb -> Z''0'  
      ELSEIF(MSTP(44).EQ.4) THEN    
        PROC(141)='f + fb -> gamma*/Z0' 
      ELSEIF(MSTP(44).EQ.5) THEN    
        PROC(141)='f + fb -> gamma*/Z''0'   
      ELSEIF(MSTP(44).EQ.6) THEN    
        PROC(141)='f + fb -> Z0/Z''0'   
      ELSEIF(MSTP(44).EQ.7) THEN    
        PROC(141)='f + fb -> gamma*/Z0/Z''0'    
      ENDIF 
    
      RETURN    
      END   
