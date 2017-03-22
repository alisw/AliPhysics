* $Id$
    
C*********************************************************************  
    
      SUBROUTINE PYOVLY_HIJING(MOVLY)  
    
C...Initializes multiplicity distribution and selects mutliplicity  
C...of overlayed events, i.e. several events occuring at the same   
C...beam crossing.  
#include "ludat1_hijing.inc"
#include "pypars_hijing.inc"
#include "pyint1_hijing.inc"
      DIMENSION WTI(0:100)  
      SAVE IMAX,WTI,WTS 
    
C...Sum of allowed cross-sections for overlayed events. 
      IF(MOVLY.EQ.1) THEN   
        VINT(131)=VINT(106) 
        IF(MSTP(132).GE.2) VINT(131)=VINT(131)+VINT(104)    
        IF(MSTP(132).GE.3) VINT(131)=VINT(131)+VINT(103)    
        IF(MSTP(132).GE.4) VINT(131)=VINT(131)+VINT(102)    
    
C...Initialize multiplicity distribution for unbiased events.   
        IF(MSTP(133).EQ.1) THEN 
          XNAVE=VINT(131)*PARP(131) 
          IF(XNAVE.GT.40.) WRITE(MSTU(11),1000) XNAVE   
          WTI(0)=EXP(-MIN(50.,XNAVE))   
          WTS=0.    
          WTN=0.    
          DO 100 I=1,100    
          WTI(I)=WTI(I-1)*XNAVE/I   
          IF(I-2.5.GT.XNAVE.AND.WTI(I).LT.1E-6) GOTO 110    
          WTS=WTS+WTI(I)    
          WTN=WTN+WTI(I)*I  
  100     IMAX=I    
  110     VINT(132)=XNAVE   
          VINT(133)=WTN/WTS 
          VINT(134)=WTS 
    
C...Initialize mutiplicity distribution for biased events.  
        ELSEIF(MSTP(133).EQ.2) THEN 
          XNAVE=VINT(131)*PARP(131) 
          IF(XNAVE.GT.40.) WRITE(MSTU(11),1000) XNAVE   
          WTI(1)=EXP(-MIN(50.,XNAVE))*XNAVE 
          WTS=WTI(1)    
          WTN=WTI(1)    
          DO 120 I=2,100    
          WTI(I)=WTI(I-1)*XNAVE/(I-1)   
          IF(I-2.5.GT.XNAVE.AND.WTI(I).LT.1E-6) GOTO 130    
          WTS=WTS+WTI(I)    
          WTN=WTN+WTI(I)*I  
  120     IMAX=I    
  130     VINT(132)=XNAVE   
          VINT(133)=WTN/WTS 
          VINT(134)=WTS 
        ENDIF   
    
C...Pick multiplicity of overlayed events.  
      ELSE  
        IF(MSTP(133).EQ.0) THEN 
          MINT(81)=MAX(1,MSTP(134)) 
        ELSE    
          WTR=WTS*RLU_HIJING(0)    
          DO 140 I=1,IMAX   
          MINT(81)=I    
          WTR=WTR-WTI(I)    
          IF(WTR.LE.0.) GOTO 150    
  140     CONTINUE  
  150     CONTINUE  
        ENDIF   
      ENDIF 
    
C...Format statement for error message. 
 1000 FORMAT(1X,'Warning: requested average number of events per bunch',    
     &'crossing too large, ',1P,E12.4)  
    
      RETURN    
      END   
