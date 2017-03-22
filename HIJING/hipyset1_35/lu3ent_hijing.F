* $Id$
    
C*********************************************************************  
    
      SUBROUTINE LU3ENT_HIJING(IP,KF1,KF2,KF3,PECM,X1,X3)  
    
C...Purpose: to store three partons or particles in their CM frame, 
C...with the first along the +z axis and the third in the (x,z) 
C...plane with x > 0.   
#include "lujets_hijing.inc"
#include "ludat1_hijing.inc"
#include "ludat2_hijing.inc"
    
C...Standard checks.    
      MSTU(28)=0    
      IF(MSTU(12).GE.1) CALL LULIST_HIJING(0)  
      IPA=MAX(1,IABS(IP))   
      IF(IPA.GT.MSTU(4)-2) CALL LUERRM_HIJING(21,  
     &'(LU3ENT_HIJING:) writing outside LUJETS_HIJING memory')    
      KC1=LUCOMP_HIJING(KF1)   
      KC2=LUCOMP_HIJING(KF2)   
      KC3=LUCOMP_HIJING(KF3)   
      IF(KC1.EQ.0.OR.KC2.EQ.0.OR.KC3.EQ.0) CALL LUERRM_HIJING(12,  
     &'(LU3ENT_HIJING:) unknown flavour code') 
    
C...Find masses. Reset K, P and V vectors.  
      PM1=0.    
      IF(MSTU(10).EQ.1) PM1=P(IPA,5)    
      IF(MSTU(10).GE.2) PM1=ULMASS_HIJING(KF1) 
      PM2=0.    
      IF(MSTU(10).EQ.1) PM2=P(IPA+1,5)  
      IF(MSTU(10).GE.2) PM2=ULMASS_HIJING(KF2) 
      PM3=0.    
      IF(MSTU(10).EQ.1) PM3=P(IPA+2,5)  
      IF(MSTU(10).GE.2) PM3=ULMASS_HIJING(KF3) 
      DO 100 I=IPA,IPA+2    
      DO 100 J=1,5  
      K(I,J)=0  
      P(I,J)=0. 
  100 V(I,J)=0. 
    
C...Check flavours. 
      KQ1=KCHG(KC1,2)*ISIGN(1,KF1)  
      KQ2=KCHG(KC2,2)*ISIGN(1,KF2)  
      KQ3=KCHG(KC3,2)*ISIGN(1,KF3)  
      IF(KQ1.EQ.0.AND.KQ2.EQ.0.AND.KQ3.EQ.0) THEN   
      ELSEIF(KQ1.NE.0.AND.KQ2.EQ.2.AND.(KQ1+KQ3.EQ.0.OR.KQ1+KQ3.EQ.4))  
     &THEN  
      ELSE  
         CALL LUERRM_HIJING(2
     $        ,'(LU3ENT_HIJING:) unphysical flavour combination')   
      ENDIF 
      K(IPA,2)=KF1  
      K(IPA+1,2)=KF2    
      K(IPA+2,2)=KF3    
    
C...Store partons/particles in K vectors for normal case.   
      IF(IP.GE.0) THEN  
        K(IPA,1)=1  
        IF(KQ1.NE.0.AND.(KQ2.NE.0.OR.KQ3.NE.0)) K(IPA,1)=2  
        K(IPA+1,1)=1    
        IF(KQ2.NE.0.AND.KQ3.NE.0) K(IPA+1,1)=2  
        K(IPA+2,1)=1    
    
C...Store partons in K vectors for parton shower evolution. 
      ELSE  
        IF(KQ1.EQ.0.OR.KQ2.EQ.0.OR.KQ3.EQ.0) CALL LUERRM_HIJING(2, 
     &        '(LU3ENT_HIJING:) requested flavours can not develop'
     $        //' parton shower')   
        K(IPA,1)=3  
        K(IPA+1,1)=3    
        K(IPA+2,1)=3    
        KCS=4   
        IF(KQ1.EQ.-1) KCS=5 
        K(IPA,KCS)=MSTU(5)*(IPA+1)  
        K(IPA,9-KCS)=MSTU(5)*(IPA+2)    
        K(IPA+1,KCS)=MSTU(5)*(IPA+2)    
        K(IPA+1,9-KCS)=MSTU(5)*IPA  
        K(IPA+2,KCS)=MSTU(5)*IPA    
        K(IPA+2,9-KCS)=MSTU(5)*(IPA+1)  
      ENDIF 
    
C...Check kinematics.   
      MKERR=0   
      IF(0.5*X1*PECM.LE.PM1.OR.0.5*(2.-X1-X3)*PECM.LE.PM2.OR.   
     &0.5*X3*PECM.LE.PM3) MKERR=1   
      PA1=SQRT(MAX(0.,(0.5*X1*PECM)**2-PM1**2)) 
      PA2=SQRT(MAX(0.,(0.5*(2.-X1-X3)*PECM)**2-PM2**2)) 
      PA3=SQRT(MAX(0.,(0.5*X3*PECM)**2-PM3**2)) 
      CTHE2=(PA3**2-PA1**2-PA2**2)/(2.*PA1*PA2) 
      CTHE3=(PA2**2-PA1**2-PA3**2)/(2.*PA1*PA3) 
      IF(ABS(CTHE2).GE.1.001.OR.ABS(CTHE3).GE.1.001) MKERR=1    
      CTHE3=MAX(-1.,MIN(1.,CTHE3))  
      IF(MKERR.NE.0) CALL LUERRM_HIJING(13,    
     &'(LU3ENT_HIJING:) unphysical kinematical variable setup')    
    
C...Store partons/particles in P vectors.   
      P(IPA,3)=PA1  
      P(IPA,4)=SQRT(PA1**2+PM1**2)  
      P(IPA,5)=PM1  
      P(IPA+2,1)=PA3*SQRT(1.-CTHE3**2)  
      P(IPA+2,3)=PA3*CTHE3  
      P(IPA+2,4)=SQRT(PA3**2+PM3**2)    
      P(IPA+2,5)=PM3    
      P(IPA+1,1)=-P(IPA+2,1)    
      P(IPA+1,3)=-P(IPA,3)-P(IPA+2,3)   
      P(IPA+1,4)=SQRT(P(IPA+1,1)**2+P(IPA+1,3)**2+PM2**2)   
      P(IPA+1,5)=PM2    
    
C...Set N. Optionally fragment/decay.   
      N=IPA+2   
      IF(IP.EQ.0) CALL LUEXEC_HIJING   
    
      RETURN    
      END   
