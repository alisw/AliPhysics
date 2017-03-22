* $Id$
    
C*********************************************************************  
    
      SUBROUTINE PYINKI_HIJING(CHFRAM,CHBEAM,CHTARG,WIN)   
    
C...Identifies the two incoming particles and sets up kinematics,   
C...including rotations and boosts to/from CM frame.    
#include "lujets_hijing.inc"
#include "ludat1_hijing.inc"
#include "pysubs_hijing.inc"
#include "pypars_hijing.inc"
#include "pyint1_hijing.inc"
      CHARACTER CHFRAM*8,CHBEAM*8,CHTARG*8,CHCOM(3)*8,CHALP(2)*26,  
     &CHIDNT(3)*8,CHTEMP*8,CHCDE(18)*8,CHINIT*76    
      DIMENSION LEN(3),KCDE(18) 
      DATA CHALP/'abcdefghijklmnopqrstuvwxyz',  
     &'ABCDEFGHIJKLMNOPQRSTUVWXYZ'/ 
      DATA CHCDE/'e-      ','e+      ','nue     ','nue~    ',   
     &'mu-     ','mu+     ','numu    ','numu~   ','tau-    ',   
     &'tau+    ','nutau   ','nutau~  ','pi+     ','pi-     ',   
     &'n       ','n~      ','p       ','p~      '/  
      DATA KCDE/11,-11,12,-12,13,-13,14,-14,15,-15,16,-16,  
     &211,-211,2112,-2112,2212,-2212/   
    
C...Convert character variables to lowercase and find their length. 
      CHCOM(1)=CHFRAM   
      CHCOM(2)=CHBEAM   
      CHCOM(3)=CHTARG   
      DO 120 I=1,3  
      LEN(I)=8  
      DO 100 LL=8,1,-1  
      IF(LEN(I).EQ.LL.AND.CHCOM(I)(LL:LL).EQ.' ') LEN(I)=LL-1   
      DO 100 LA=1,26    
  100 IF(CHCOM(I)(LL:LL).EQ.CHALP(2)(LA:LA)) CHCOM(I)(LL:LL)=   
     &CHALP(1)(LA:LA)   
      CHIDNT(I)=CHCOM(I)    
      DO 110 LL=1,6 
      IF(CHIDNT(I)(LL:LL+2).EQ.'bar') THEN  
        CHTEMP=CHIDNT(I)    
        CHIDNT(I)=CHTEMP(1:LL-1)//'~'//CHTEMP(LL+3:8)//'  ' 
      ENDIF 
  110 CONTINUE  
      DO 120 LL=1,8 
      IF(CHIDNT(I)(LL:LL).EQ.'_') THEN  
        CHTEMP=CHIDNT(I)    
        CHIDNT(I)=CHTEMP(1:LL-1)//CHTEMP(LL+1:8)//' '   
      ENDIF 
  120 CONTINUE  
    
C...Set initial state. Error for unknown codes. Reset variables.    
      N=2   
      DO 140 I=1,2  
      K(I,2)=0  
      DO 130 J=1,18 
  130 IF(CHIDNT(I+1).EQ.CHCDE(J)) K(I,2)=KCDE(J)    
      P(I,5)=ULMASS_HIJING(K(I,2)) 
      MINT(40+I)=1  
      IF(IABS(K(I,2)).GT.100) MINT(40+I)=2  
      DO 140 J=1,5  
  140 V(I,J)=0. 
      IF(K(1,2).EQ.0) WRITE(MSTU(11),1000) CHBEAM(1:LEN(2)) 
      IF(K(2,2).EQ.0) WRITE(MSTU(11),1100) CHTARG(1:LEN(3)) 
      IF(K(1,2).EQ.0.OR.K(2,2).EQ.0) STOP   
      DO 150 J=6,10 
  150 VINT(J)=0.    
      CHINIT=' '    
    
C...Set up kinematics for events defined in CM frame.   
      IF(CHCOM(1)(1:2).EQ.'cm') THEN    
        IF(CHCOM(2)(1:1).NE.'e') THEN   
          LOFFS=(34-(LEN(2)+LEN(3)))/2  
          CHINIT(LOFFS+1:76)='PYTHIA will be initialized for a '//  
     &    CHCOM(2)(1:LEN(2))//'-'//CHCOM(3)(1:LEN(3))//' collider'//' ' 
        ELSE    
          LOFFS=(33-(LEN(2)+LEN(3)))/2  
          CHINIT(LOFFS+1:76)='PYTHIA will be initialized for an '// 
     &    CHCOM(2)(1:LEN(2))//'-'//CHCOM(3)(1:LEN(3))//' collider'//' ' 
        ENDIF   
C        WRITE(MSTU(11),1200) CHINIT 
C        WRITE(MSTU(11),1300) WIN    
        S=WIN**2    
        P(1,1)=0.   
        P(1,2)=0.   
        P(2,1)=0.   
        P(2,2)=0.   
        P(1,3)=SQRT(((S-P(1,5)**2-P(2,5)**2)**2-(2.*P(1,5)*P(2,5))**2)/ 
     &  (4.*S)) 
        P(2,3)=-P(1,3)  
        P(1,4)=SQRT(P(1,3)**2+P(1,5)**2)    
        P(2,4)=SQRT(P(2,3)**2+P(2,5)**2)    
    
C...Set up kinematics for fixed target events.  
      ELSEIF(CHCOM(1)(1:3).EQ.'fix') THEN   
        LOFFS=(29-(LEN(2)+LEN(3)))/2    
        CHINIT(LOFFS+1:76)='PYTHIA will be initialized for '//  
     &  CHCOM(2)(1:LEN(2))//' on '//CHCOM(3)(1:LEN(3))//    
     &  ' fixed target'//' '    
C        WRITE(MSTU(11),1200) CHINIT 
C        WRITE(MSTU(11),1400) WIN    
        P(1,1)=0.   
        P(1,2)=0.   
        P(2,1)=0.   
        P(2,2)=0.   
        P(1,3)=WIN  
        P(1,4)=SQRT(P(1,3)**2+P(1,5)**2)    
        P(2,3)=0.   
        P(2,4)=P(2,5)   
        S=P(1,5)**2+P(2,5)**2+2.*P(2,4)*P(1,4)  
        VINT(10)=P(1,3)/(P(1,4)+P(2,4)) 
        CALL LUROBO_HIJING(0.,0.,0.,0.,-VINT(10))  
C        WRITE(MSTU(11),1500) SQRT(S)    
    
C...Set up kinematics for events in user-defined frame. 
      ELSEIF(CHCOM(1)(1:3).EQ.'use') THEN   
        LOFFS=(13-(LEN(1)+LEN(2)))/2    
        CHINIT(LOFFS+1:76)='PYTHIA will be initialized for '//  
     &  CHCOM(2)(1:LEN(2))//' on '//CHCOM(3)(1:LEN(3))//    
     &  'user-specified configuration'//' ' 
C        WRITE(MSTU(11),1200) CHINIT 
C        WRITE(MSTU(11),1600)    
C        WRITE(MSTU(11),1700) CHCOM(2),P(1,1),P(1,2),P(1,3)  
C        WRITE(MSTU(11),1700) CHCOM(3),P(2,1),P(2,2),P(2,3)  
        P(1,4)=SQRT(P(1,1)**2+P(1,2)**2+P(1,3)**2+P(1,5)**2)    
        P(2,4)=SQRT(P(2,1)**2+P(2,2)**2+P(2,3)**2+P(2,5)**2)    
        DO 160 J=1,3    
  160   VINT(7+J)=(DBLE(P(1,J))+DBLE(P(2,J)))/DBLE(P(1,4)+P(2,4))   
        CALL LUROBO_HIJING(0.,0.,-VINT(8),-VINT(9),-VINT(10))  
        VINT(7)=ULANGL_HIJING(P(1,1),P(1,2))   
        CALL LUROBO_HIJING(0.,-VINT(7),0.,0.,0.)   
        VINT(6)=ULANGL_HIJING(P(1,3),P(1,1))   
        CALL LUROBO_HIJING(-VINT(6),0.,0.,0.,0.)   
        S=P(1,5)**2+P(2,5)**2+2.*(P(1,4)*P(2,4)-P(1,3)*P(2,3))  
C        WRITE(MSTU(11),1500) SQRT(S)    
    
C...Unknown frame. Error for too low CM energy. 
      ELSE  
        WRITE(MSTU(11),1800) CHFRAM(1:LEN(1))   
        STOP    
      ENDIF 
      IF(S.LT.PARP(2)**2) THEN  
        WRITE(MSTU(11),1900) SQRT(S)    
        STOP    
      ENDIF 
    
C...Save information on incoming particles. 
      MINT(11)=K(1,2)   
      MINT(12)=K(2,2)   
      MINT(43)=2*MINT(41)+MINT(42)-2    
      VINT(1)=SQRT(S)   
      VINT(2)=S 
      VINT(3)=P(1,5)    
      VINT(4)=P(2,5)    
      VINT(5)=P(1,3)    
    
C...Store constants to be used in generation.   
      IF(MSTP(82).LE.1) VINT(149)=4.*PARP(81)**2/S  
      IF(MSTP(82).GE.2) VINT(149)=4.*PARP(82)**2/S  
    
C...Formats for initialization and error information.   
 1000 FORMAT(1X,'Error: unrecognized beam particle ''',A,'''.'/ 
     &1X,'Execution stopped!')  
 1100 FORMAT(1X,'Error: unrecognized target particle ''',A,'''.'/   
     &1X,'Execution stopped!')  
 1200 FORMAT(/1X,78('=')/1X,'I',76X,'I'/1X,'I',A76,'I') 
 1300 FORMAT(1X,'I',18X,'at',1X,F10.3,1X,'GeV center-of-mass energy',   
     &19X,'I'/1X,'I',76X,'I'/1X,78('='))    
 1400 FORMAT(1X,'I',22X,'at',1X,F10.3,1X,'GeV/c lab-momentum',22X,'I')  
 1500 FORMAT(1X,'I',76X,'I'/1X,'I',11X,'corresponding to',1X,F10.3,1X,  
     &'GeV center-of-mass energy',12X,'I'/1X,'I',76X,'I'/1X,78('='))    
 1600 FORMAT(1X,'I',76X,'I'/1X,'I',24X,'px (GeV/c)',3X,'py (GeV/c)',3X, 
     &'pz (GeV/c)',16X,'I') 
 1700 FORMAT(1X,'I',15X,A8,3(2X,F10.3,1X),15X,'I')  
 1800 FORMAT(1X,'Error: unrecognized coordinate frame ''',A,'''.'/  
     &1X,'Execution stopped!')  
 1900 FORMAT(1X,'Error: too low CM energy,',F8.3,' GeV for event ', 
     &'generation.'/1X,'Execution stopped!')    
    
      RETURN    
      END   
