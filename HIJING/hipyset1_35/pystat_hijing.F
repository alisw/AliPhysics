* $Id$
    
C***********************************************************************    
    
      SUBROUTINE PYSTAT_HIJING(MSTAT)  
    
C...Prints out information about cross-sections, decay widths, branching    
C...ratios, kinematical limits, status codes and parameter values.  
#include "ludat1_hijing.inc"
#include "ludat2_hijing.inc"
#include "ludat3_hijing.inc"
#include "pysubs_hijing.inc"
#include "pypars_hijing.inc"
#include "pyint1_hijing.inc"
#include "pyint4_hijing.inc"
#include "pyint5_hijing.inc"
#include "pyint6_hijing.inc"
      CHARACTER CHAU*16,CHPA(-40:40)*12,CHIN(2)*12, 
     &STATE(-1:5)*4,CHKIN(21)*18    
      DATA STATE/'----','off ','on  ','on/+','on/-','on/1','on/2'/, 
     &CHKIN/' m_hard (GeV/c^2) ',' p_T_hard (GeV/c) ',  
     &'m_finite (GeV/c^2)','   y*_subsystem   ','     y*_large     ',   
     &'     y*_small     ','    eta*_large    ','    eta*_small    ',   
     &'cos(theta*)_large ','cos(theta*)_small ','       x_1        ',   
     &'       x_2        ','       x_F        ',' cos(theta_hard)  ',   
     &'m''_hard (GeV/c^2) ','       tau        ','        y*        ',  
     &'cos(theta_hard^-) ','cos(theta_hard^+) ','      x_T^2       ',   
     &'       tau''       '/    
    
C...Cross-sections. 
      IF(MSTAT.LE.1) THEN   
        WRITE(MSTU(11),1000)    
        WRITE(MSTU(11),1100)    
        WRITE(MSTU(11),1200) 0,PROC(0),NGEN(0,3),NGEN(0,1),XSEC(0,3)    
        DO 100 I=1,200  
        IF(MSUB(I).NE.1) GOTO 100   
        WRITE(MSTU(11),1200) I,PROC(I),NGEN(I,3),NGEN(I,1),XSEC(I,3)    
  100   CONTINUE    
        WRITE(MSTU(11),1300) 1.-FLOAT(NGEN(0,3))/   
     &  MAX(1.,FLOAT(NGEN(0,2)))    
    
C...Decay widths and branching ratios.  
      ELSEIF(MSTAT.EQ.2) THEN   
        DO 110 KF=-40,40    
        CALL LUNAME_HIJING(KF,CHAU)    
  110   CHPA(KF)=CHAU(1:12) 
        WRITE(MSTU(11),1400)    
        WRITE(MSTU(11),1500)    
C...Off-shell branchings.   
        DO 130 I=1,17   
        KC=I    
        IF(I.GE.9) KC=I+2   
        IF(I.EQ.17) KC=21   
        WRITE(MSTU(11),1600) CHPA(KC),0.,0.,STATE(MDCY(KC,1)),0.    
        DO 120 J=1,MDCY(KC,3)   
        IDC=J+MDCY(KC,2)-1  
  120   IF(MDME(IDC,2).EQ.102) WRITE(MSTU(11),1700) CHPA(KFDP(IDC,1)),  
     &  CHPA(KFDP(IDC,2)),0.,0.,STATE(MDME(IDC,1)),0.   
  130   CONTINUE    
C...On-shell decays.    
        DO 150 I=1,6    
        KC=I+22 
        IF(I.EQ.4) KC=32    
        IF(I.EQ.5) KC=37    
        IF(I.EQ.6) KC=40    
        IF(WIDE(KC,0).GT.0.) THEN   
          WRITE(MSTU(11),1600) CHPA(KC),WIDP(KC,0),1.,  
     &    STATE(MDCY(KC,1)),1.  
          DO 140 J=1,MDCY(KC,3) 
          IDC=J+MDCY(KC,2)-1    
  140     WRITE(MSTU(11),1700) CHPA(KFDP(IDC,1)),CHPA(KFDP(IDC,2)), 
     &    WIDP(KC,J),WIDP(KC,J)/WIDP(KC,0),STATE(MDME(IDC,1)),  
     &    WIDE(KC,J)/WIDE(KC,0) 
        ELSE    
          WRITE(MSTU(11),1600) CHPA(KC),WIDP(KC,0),1.,  
     &    STATE(MDCY(KC,1)),0.  
        ENDIF   
  150   CONTINUE    
        WRITE(MSTU(11),1800)    
    
C...Allowed incoming partons/particles at hard interaction. 
      ELSEIF(MSTAT.EQ.3) THEN   
        WRITE(MSTU(11),1900)    
        CALL LUNAME_HIJING(MINT(11),CHAU)  
        CHIN(1)=CHAU(1:12)  
        CALL LUNAME_HIJING(MINT(12),CHAU)  
        CHIN(2)=CHAU(1:12)  
        WRITE(MSTU(11),2000) CHIN(1),CHIN(2)    
        DO 160 KF=-40,40    
        CALL LUNAME_HIJING(KF,CHAU)    
  160   CHPA(KF)=CHAU(1:12) 
        IF(MINT(43).EQ.1) THEN  
          WRITE(MSTU(11),2100) CHPA(MINT(11)),STATE(KFIN(1,MINT(11))),  
     &    CHPA(MINT(12)),STATE(KFIN(2,MINT(12)))    
        ELSEIF(MINT(43).EQ.2) THEN  
          WRITE(MSTU(11),2100) CHPA(MINT(11)),STATE(KFIN(1,MINT(11))),  
     &    CHPA(-MSTP(54)),STATE(KFIN(2,-MSTP(54)))  
          DO 170 I=-MSTP(54)+1,-1   
  170     WRITE(MSTU(11),2200) CHPA(I),STATE(KFIN(2,I)) 
          DO 180 I=1,MSTP(54)   
  180     WRITE(MSTU(11),2200) CHPA(I),STATE(KFIN(2,I)) 
          WRITE(MSTU(11),2200) CHPA(21),STATE(KFIN(2,21))   
        ELSEIF(MINT(43).EQ.3) THEN  
          WRITE(MSTU(11),2100) CHPA(-MSTP(54)),STATE(KFIN(1,-MSTP(54))),    
     &    CHPA(MINT(12)),STATE(KFIN(2,MINT(12)))    
          DO 190 I=-MSTP(54)+1,-1   
  190     WRITE(MSTU(11),2300) CHPA(I),STATE(KFIN(1,I)) 
          DO 200 I=1,MSTP(54)   
  200     WRITE(MSTU(11),2300) CHPA(I),STATE(KFIN(1,I)) 
          WRITE(MSTU(11),2300) CHPA(21),STATE(KFIN(1,21))   
        ELSEIF(MINT(43).EQ.4) THEN  
          DO 210 I=-MSTP(54),-1 
  210     WRITE(MSTU(11),2100) CHPA(I),STATE(KFIN(1,I)),CHPA(I),    
     &    STATE(KFIN(2,I))  
          DO 220 I=1,MSTP(54)   
  220     WRITE(MSTU(11),2100) CHPA(I),STATE(KFIN(1,I)),CHPA(I),    
     &    STATE(KFIN(2,I))  
          WRITE(MSTU(11),2100) CHPA(21),STATE(KFIN(1,21)),CHPA(21), 
     &    STATE(KFIN(2,21)) 
        ENDIF   
        WRITE(MSTU(11),2400)    
    
C...User-defined and derived limits on kinematical variables.   
      ELSEIF(MSTAT.EQ.4) THEN   
        WRITE(MSTU(11),2500)    
        WRITE(MSTU(11),2600)    
        SHRMAX=CKIN(2)  
        IF(SHRMAX.LT.0.) SHRMAX=VINT(1) 
        WRITE(MSTU(11),2700) CKIN(1),CHKIN(1),SHRMAX    
        PTHMIN=MAX(CKIN(3),CKIN(5)) 
        PTHMAX=CKIN(4)  
        IF(PTHMAX.LT.0.) PTHMAX=0.5*SHRMAX  
        WRITE(MSTU(11),2800) CKIN(3),PTHMIN,CHKIN(2),PTHMAX 
        WRITE(MSTU(11),2900) CHKIN(3),CKIN(6)   
        DO 230 I=4,14   
  230   WRITE(MSTU(11),2700) CKIN(2*I-1),CHKIN(I),CKIN(2*I) 
        SPRMAX=CKIN(32) 
        IF(SPRMAX.LT.0.) SPRMAX=VINT(1) 
        WRITE(MSTU(11),2700) CKIN(31),CHKIN(13),SPRMAX  
        WRITE(MSTU(11),3000)    
        WRITE(MSTU(11),3100)    
        WRITE(MSTU(11),2600)    
        DO 240 I=16,21  
  240   WRITE(MSTU(11),2700) VINT(I-5),CHKIN(I),VINT(I+15)  
        WRITE(MSTU(11),3000)    
    
C...Status codes and parameter values.  
      ELSEIF(MSTAT.EQ.5) THEN   
        WRITE(MSTU(11),3200)    
        WRITE(MSTU(11),3300)    
        DO 250 I=1,100  
  250   WRITE(MSTU(11),3400) I,MSTP(I),PARP(I),100+I,MSTP(100+I),   
     &  PARP(100+I) 
      ENDIF 
    
C...Formats for printouts.  
 1000 FORMAT('1',9('*'),1X,'PYSTAT_HIJING:  Statistics on Number of ', 
     &'Events and Cross-sections',1X,9('*'))    
 1100 FORMAT(/1X,78('=')/1X,'I',34X,'I',28X,'I',12X,'I'/1X,'I',12X, 
     &'Subprocess',12X,'I',6X,'Number of points',6X,'I',4X,'Sigma',3X,  
     &'I'/1X,'I',34X,'I',28X,'I',12X,'I'/1X,'I',34('-'),'I',28('-'),    
     &'I',4X,'(mb)',4X,'I'/1X,'I',34X,'I',28X,'I',12X,'I'/1X,'I',1X,    
     &'N:o',1X,'Type',25X,'I',4X,'Generated',9X,'Tried',1X,'I',12X, 
     &'I'/1X,'I',34X,'I',28X,'I',12X,'I'/1X,78('=')/1X,'I',34X,'I',28X, 
     &'I',12X,'I')  
 1200 FORMAT(1X,'I',1X,I3,1X,A28,1X,'I',1X,I12,1X,I13,1X,'I',1X,1P, 
     &E10.3,1X,'I') 
 1300 FORMAT(1X,'I',34X,'I',28X,'I',12X,'I'/1X,78('=')//    
     &1X,'********* Fraction of events that fail fragmentation ',   
     &'cuts =',1X,F8.5,' *********'/)   
 1400 FORMAT('1',17('*'),1X
     $     ,'PYSTAT_HIJING:  Decay Widths and Branching ','Ratios',1X
     $     ,17('*'))  
 1500 FORMAT(/1X,78('=')/1X,'I',29X,'I',13X,'I',12X,'I',6X,'I',12X,'I'/ 
     &1X,'I',1X,'Branching/Decay Channel',5X,'I',1X,'Width (GeV)',1X,   
     &'I',7X,'B.R.',1X,'I',1X,'Stat',1X,'I',2X,'Eff. B.R.',1X,'I'/1X,   
     &'I',29X,'I',13X,'I',12X,'I',6X,'I',12X,'I'/1X,78('='))    
 1600 FORMAT(1X,'I',29X,'I',13X,'I',12X,'I',6X,'I',12X,'I'/1X,'I',1X,   
     &A12,1X,'->',13X,'I',2X,1P,E10.3,0P,1X,'I',1X,1P,E10.3,0P,1X,'I',  
     &1X,A4,1X,'I',1X,1P,E10.3,0P,1X,'I')   
 1700 FORMAT(1X,'I',1X,A12,1X,'+',1X,A12,1X,'I',2X,1P,E10.3,0P,1X,'I',  
     &1X,1P,E10.3,0P,1X,'I',1X,A4,1X,'I',1X,1P,E10.3,0P,1X,'I') 
 1800 FORMAT(1X,'I',29X,'I',13X,'I',12X,'I',6X,'I',12X,'I'/1X,78('='))  
 1900 FORMAT('1',7('*'),1X,'PYSTAT_HIJING: Allowed Incoming Partons/', 
     &'Particles at Hard Interaction',1X,7('*'))    
 2000 FORMAT(/1X,78('=')/1X,'I',38X,'I',37X,'I'/1X,'I',1X,  
     &'Beam particle:',1X,A,10X,'I',1X,'Target particle:',1X,A,7X,  
     &'I'/1X,'I',38X,'I',37X,'I'/1X,'I',1X,'Content',9X,'State',16X,    
     &'I',1X,'Content',9X,'State',15X,'I'/1X,'I',38X,'I',37X,'I'/1X,    
     &78('=')/1X,'I',38X,'I',37X,'I')   
 2100 FORMAT(1X,'I',1X,A,5X,A,16X,'I',1X,A,5X,A,15X,'I')    
 2200 FORMAT(1X,'I',38X,'I',1X,A,5X,A,15X,'I')  
 2300 FORMAT(1X,'I',1X,A,5X,A,16X,'I',37X,'I')  
 2400 FORMAT(1X,'I',38X,'I',37X,'I'/1X,78('=')) 
 2500 FORMAT('1',12('*'),1X,'PYSTAT_HIJING: User-Defined Limits on ',  
     &'Kinematical Variables',1X,12('*'))   
 2600 FORMAT(/1X,78('=')/1X,'I',76X,'I')    
 2700 FORMAT(1X,'I',16X,1P,E10.3,0P,1X,'<',1X,A,1X,'<',1X,1P,E10.3,0P,  
     &16X,'I')  
 2800 FORMAT(1X,'I',3X,1P,E10.3,0P,1X,'(',1P,E10.3,0P,')',1X,'<',1X,A,  
     &1X,'<',1X,1P,E10.3,0P,16X,'I')    
 2900 FORMAT(1X,'I',29X,A,1X,'=',1X,1P,E10.3,0P,16X,'I')    
 3000 FORMAT(1X,'I',76X,'I'/1X,78('=')) 
 3100 FORMAT(////1X,5('*'),1X
     $     ,'PYSTAT_HIJING: Derived Limits on Kinematical '
     $     ,'Variables Used in Generation',1X,5('*')) 
 3200 FORMAT('1',12('*'),1X
     $     ,'PYSTAT_HIJING: Summary of Status Codes and '
     $     ,'Parameter Values',1X,12('*'))    
 3300 FORMAT(/3X,'I',4X,'MSTP(I)',9X,'PARP(I)',20X,'I',4X,'MSTP(I)',9X, 
     &'PARP(I)'/)   
 3400 FORMAT(1X,I3,5X,I6,6X,1P,E10.3,0P,18X,I3,5X,I6,6X,1P,E10.3)   
    
      RETURN    
      END   
