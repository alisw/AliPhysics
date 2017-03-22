* $Id$
C
C
C
C     Modified for HIJING program
c
c    modification July 22, 1997  In pyremnn put an upper limit
c     on the total pt kick the parton can accumulate via multiple
C     scattering. Set the upper limit to be the sqrt(s)/2,
c     this is fix cronin bug for Pb+Pb events at SPS energy.
c
C
C Last modification Oct. 1993 to comply with non-vax
C machines' compiler 
C
C
      SUBROUTINE LU1ENT_HIJING(IP,KF,PE,THE,PHI)   
    
C...Purpose: to store one parton/particle in commonblock LUJETS.    
#include "lujets_hijing.inc"
#include "ludat1_hijing.inc"
#include "ludat2_hijing.inc"
    
C...Standard checks.    
      MSTU(28)=0    
      IF(MSTU(12).GE.1) CALL LULIST_HIJING(0)  
      IPA=MAX(1,IABS(IP))   
      IF(IPA.GT.MSTU(4)) CALL LUERRM_HIJING(21,    
     &'(LU1ENT_HIJING:) writing outside LUJETS_HIJING memory')    
      KC=LUCOMP_HIJING(KF) 
      IF(KC.EQ.0) CALL LUERRM_HIJING(12
     $     ,'(LU1ENT_HIJING:) unknown flavour code')  
    
C...Find mass. Reset K, P and V vectors.    
      PM=0. 
      IF(MSTU(10).EQ.1) PM=P(IPA,5) 
      IF(MSTU(10).GE.2) PM=ULMASS_HIJING(KF)   
      DO 100 J=1,5  
      K(IPA,J)=0    
      P(IPA,J)=0.   
  100 V(IPA,J)=0.   
    
C...Store parton/particle in K and P vectors.   
      K(IPA,1)=1    
      IF(IP.LT.0) K(IPA,1)=2    
      K(IPA,2)=KF   
      P(IPA,5)=PM   
      P(IPA,4)=MAX(PE,PM)   
      PA=SQRT(P(IPA,4)**2-P(IPA,5)**2)  
      P(IPA,1)=PA*SIN(THE)*COS(PHI) 
      P(IPA,2)=PA*SIN(THE)*SIN(PHI) 
      P(IPA,3)=PA*COS(THE)  
    
C...Set N. Optionally fragment/decay.   
      N=IPA 
      IF(IP.EQ.0) CALL LUEXEC_HIJING   
    
      RETURN    
      END   
