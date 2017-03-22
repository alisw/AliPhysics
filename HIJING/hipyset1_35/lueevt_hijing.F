* $Id$
    
C*********************************************************************  
    
      SUBROUTINE LUEEVT_HIJING(KFL,ECM)    
    
C...Purpose: to handle the generation of an e+e- annihilation jet event.    
      IMPLICIT DOUBLE PRECISION(D)  
#include "lujets_hijing.inc"
#include "ludat1_hijing.inc"
#include "ludat2_hijing.inc"
    
C...Check input parameters. 
      IF(MSTU(12).GE.1) CALL LULIST_HIJING(0)  
      IF(KFL.LT.0.OR.KFL.GT.8) THEN 
         CALL LUERRM_HIJING(16
     $        ,'(LUEEVT_HIJING:) called with unknown flavour code')    
        IF(MSTU(21).GE.1) RETURN    
      ENDIF 
      IF(KFL.LE.5) ECMMIN=PARJ(127)+2.02*PARF(100+MAX(1,KFL))   
      IF(KFL.GE.6) ECMMIN=PARJ(127)+2.02*PMAS(KFL,1)    
      IF(ECM.LT.ECMMIN) THEN    
         CALL LUERRM_HIJING(16
     $        ,'(LUEEVT_HIJING:) called with too small CM energy') 
        IF(MSTU(21).GE.1) RETURN    
      ENDIF 
    
C...Check consistency of MSTJ options set.  
      IF(MSTJ(109).EQ.2.AND.MSTJ(110).NE.1) THEN    
        CALL LUERRM_HIJING(6,  
     &  '(LUEEVT_HIJING:) MSTJ(109) value requires MSTJ(110) = 1') 
        MSTJ(110)=1 
      ENDIF 
      IF(MSTJ(109).EQ.2.AND.MSTJ(111).NE.0) THEN    
        CALL LUERRM_HIJING(6,  
     &  '(LUEEVT_HIJING:) MSTJ(109) value requires MSTJ(111) = 0') 
        MSTJ(111)=0 
      ENDIF 
    
C...Initialize alpha_strong and total cross-section.    
      MSTU(111)=MSTJ(108)   
      IF(MSTJ(108).EQ.2.AND.(MSTJ(101).EQ.0.OR.MSTJ(101).EQ.1)) 
     &MSTU(111)=1   
      PARU(112)=PARJ(121)   
      IF(MSTU(111).EQ.2) PARU(112)=PARJ(122)    
      IF(MSTJ(116).GT.0.AND.(MSTJ(116).GE.2.OR.ABS(ECM-PARJ(151)).GE.   
     &     PARJ(139).OR.10*MSTJ(102)+KFL.NE.MSTJ(119))) CALL
     $     LUXTOT_HIJING(KFL,ECM,XTOT) 
      IF(MSTJ(116).GE.3) MSTJ(116)=1    
    
C...Add initial e+e- to event record (documentation only).  
      NTRY=0    
  100 NTRY=NTRY+1   
      IF(NTRY.GT.100) THEN  
         CALL LUERRM_HIJING(14
     $        ,'(LUEEVT_HIJING:) caught in an infinite loop')  
        RETURN  
      ENDIF 
      NC=0  
      IF(MSTJ(115).GE.2) THEN   
        NC=NC+2 
        CALL LU1ENT_HIJING(NC-1,11,0.5*ECM,0.,0.)  
        K(NC-1,1)=21    
        CALL LU1ENT_HIJING(NC,-11,0.5*ECM,PARU(1),0.)  
        K(NC,1)=21  
      ENDIF 
    
C...Radiative photon (in initial state).    
      MK=0  
      ECMC=ECM  
      IF(MSTJ(107).GE.1.AND.MSTJ(116).GE.1) CALL LURADK_HIJING(ECM,MK
     $     ,PAK,THEK,PHIK,ALPK)   
      IF(MK.EQ.1) ECMC=SQRT(ECM*(ECM-2.*PAK))   
      IF(MSTJ(115).GE.1.AND.MK.EQ.1) THEN   
        NC=NC+1 
        CALL LU1ENT_HIJING(NC,22,PAK,THEK,PHIK)    
        K(NC,3)=MIN(MSTJ(115)/2,1)  
      ENDIF 
    
C...Virtual exchange boson (gamma or Z0).   
      IF(MSTJ(115).GE.3) THEN   
        NC=NC+1 
        KF=22   
        IF(MSTJ(102).EQ.2) KF=23    
        MSTU10=MSTU(10) 
        MSTU(10)=1  
        P(NC,5)=ECMC    
        CALL LU1ENT_HIJING(NC,KF,ECMC,0.,0.)   
        K(NC,1)=21  
        K(NC,3)=1   
        MSTU(10)=MSTU10 
      ENDIF 
    
C...Choice of flavour and jet configuration.    
      CALL LUXKFL_HIJING(KFL,ECM,ECMC,KFLC)    
      IF(KFLC.EQ.0) GOTO 100    
      CALL LUXJET_HIJING(ECMC,NJET,CUT)    
      KFLN=21   
      IF(NJET.EQ.4) CALL LUX4JT_HIJING(NJET,CUT,KFLC,ECMC,KFLN,X1,X2,X4,   
     &X12,X14)  
      IF(NJET.EQ.3) CALL LUX3JT_HIJING(NJET,CUT,KFLC,ECMC,X1,X3)   
      IF(NJET.EQ.2) MSTJ(120)=1 
    
C...Fill jet configuration and origin.  
      IF(NJET.EQ.2.AND.MSTJ(101).NE.5) CALL LU2ENT_HIJING(NC+1,KFLC,
     $     -KFLC,ECMC)    
      IF(NJET.EQ.2.AND.MSTJ(101).EQ.5) CALL LU2ENT_HIJING(-(NC+1),KFLC,
     $     -KFLC,ECMC) 
      IF(NJET.EQ.3) CALL LU3ENT_HIJING(NC+1,KFLC,21,-KFLC,ECMC,X1,X3)  
      IF(NJET.EQ.4.AND.KFLN.EQ.21) CALL LU4ENT_HIJING(NC+1,KFLC,KFLN
     $     ,KFLN,-KFLC,ECMC,X1,X2,X4,X12,X14)  
      IF(NJET.EQ.4.AND.KFLN.NE.21) CALL LU4ENT_HIJING(NC+1,KFLC,-KFLN
     $     ,KFLN,-KFLC,ECMC,X1,X2,X4,X12,X14)  
      DO 110 IP=NC+1,N  
  110 K(IP,3)=K(IP,3)+MIN(MSTJ(115)/2,1)+(MSTJ(115)/3)*(NC-1)   
    
C...Angular orientation according to matrix element.    
      IF(MSTJ(106).EQ.1) THEN   
        CALL LUXDIF_HIJING(NC,NJET,KFLC,ECMC,CHI,THE,PHI)  
        CALL LUDBRB_HIJING(NC+1,N,0.,CHI,0D0,0D0,0D0)  
        CALL LUDBRB_HIJING(NC+1,N,THE,PHI,0D0,0D0,0D0) 
      ENDIF 
    
C...Rotation and boost from radiative photon.   
      IF(MK.EQ.1) THEN  
        DBEK=-PAK/(ECM-PAK) 
        NMIN=NC+1-MSTJ(115)/3   
        CALL LUDBRB_HIJING(NMIN,N,0.,-PHIK,0D0,0D0,0D0)    
        CALL LUDBRB_HIJING(NMIN,N,ALPK,0.,DBEK*SIN(THEK),0D0,DBEK
     $       *COS(THEK))   
        CALL LUDBRB_HIJING(NMIN,N,0.,PHIK,0D0,0D0,0D0) 
      ENDIF 
    
C...Generate parton shower. Rearrange along strings and check.  
      IF(MSTJ(101).EQ.5) THEN   
        CALL LUSHOW_HIJING(N-1,N,ECMC) 
        MSTJ14=MSTJ(14) 
        IF(MSTJ(105).EQ.-1) MSTJ(14)=0  
        IF(MSTJ(105).GE.0) MSTU(28)=0   
        CALL LUPREP_HIJING(0)  
        MSTJ(14)=MSTJ14 
        IF(MSTJ(105).GE.0.AND.MSTU(28).NE.0) GOTO 100   
      ENDIF 
    
C...Fragmentation/decay generation. Information for LUTABU_HIJING. 
      IF(MSTJ(105).EQ.1) CALL LUEXEC_HIJING    
      MSTU(161)=KFLC    
      MSTU(162)=-KFLC   
    
      RETURN    
      END   
