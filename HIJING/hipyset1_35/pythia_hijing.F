* $Id$
    
C*********************************************************************  
    
      SUBROUTINE PYTHIA_HIJING 
    
C...Administers the generation of a high-pt event via calls to a number 
C...of subroutines; also computes cross-sections.   
#include "lujets_hijing.inc"
#include "ludat1_hijing.inc"
#include "ludat2_hijing.inc"
#include "pysubs_hijing.inc"
#include "pypars_hijing.inc"
#include "pyint1_hijing.inc"
#include "pyint2_hijing.inc"
#include "pyint5_hijing.inc"
    
C...Loop over desired number of overlayed events (normally 1).  
      MINT(7)=0 
      MINT(8)=0 
      NOVL=1    
      IF(MSTP(131).NE.0) CALL PYOVLY_HIJING(2) 
      IF(MSTP(131).NE.0) NOVL=MINT(81)  
      MINT(83)=0    
      MINT(84)=MSTP(126)    
      MSTU(70)=0    
      DO 190 IOVL=1,NOVL    
      IF(MINT(84)+100.GE.MSTU(4)) THEN  
        CALL LUERRM_HIJING(11, 
     &        '(PYTHIA_HIJING:) no more space in LUJETS_HIJING '/
     $        /'for overlayed events')   
        IF(MSTU(21).GE.1) GOTO 200  
      ENDIF 
      MINT(82)=IOVL 
    
C...Generate variables of hard scattering.  
  100 CONTINUE  
      IF(IOVL.EQ.1) NGEN(0,2)=NGEN(0,2)+1   
      MINT(31)=0    
      MINT(51)=0    
      CALL PYRAND_HIJING   
      ISUB=MINT(1)  
      IF(IOVL.EQ.1) THEN    
        NGEN(ISUB,2)=NGEN(ISUB,2)+1 
    
C...Store information on hard interaction.  
        DO 110 J=1,200  
        MSTI(J)=0   
  110   PARI(J)=0.  
        MSTI(1)=MINT(1) 
        MSTI(2)=MINT(2) 
        MSTI(11)=MINT(11)   
        MSTI(12)=MINT(12)   
        MSTI(15)=MINT(15)   
        MSTI(16)=MINT(16)   
        MSTI(17)=MINT(17)   
        MSTI(18)=MINT(18)   
        PARI(11)=VINT(1)    
        PARI(12)=VINT(2)    
        IF(ISUB.NE.95) THEN 
          DO 120 J=13,22    
  120     PARI(J)=VINT(30+J)    
          PARI(33)=VINT(41) 
          PARI(34)=VINT(42) 
          PARI(35)=PARI(33)-PARI(34)    
          PARI(36)=VINT(21) 
          PARI(37)=VINT(22) 
          PARI(38)=VINT(26) 
          PARI(41)=VINT(23) 
        ENDIF   
      ENDIF 
    
      IF(MSTP(111).EQ.-1) GOTO 160  
      IF(ISUB.LE.90.OR.ISUB.GE.95) THEN 
C...Hard scattering (including low-pT): 
C...reconstruct kinematics and colour flow of hard scattering.  
        CALL PYSCAT_HIJING 
        IF(MINT(51).EQ.1) GOTO 100  
    
C...Showering of initial state partons (optional).  
        IPU1=MINT(84)+1 
        IPU2=MINT(84)+2 
        IF(MSTP(61).GE.1.AND.MINT(43).NE.1.AND.ISUB.NE.95)  
     &  CALL PYSSPA_HIJING(IPU1,IPU2)  
        NSAV1=N 
    
C...Multiple interactions.  
        IF(MSTP(81).GE.1.AND.MINT(43).EQ.4.AND.ISUB.NE.95)  
     &  CALL PYMULT_HIJING(6)  
        MINT(1)=ISUB    
        NSAV2=N 
    
C...Hadron remnants and primordial kT.  
        CALL PYREMN_HIJING(IPU1,IPU2)  
        IF(MINT(51).EQ.1) GOTO 100  
        NSAV3=N 
    
C...Showering of final state partons (optional).    
        IPU3=MINT(84)+3 
        IPU4=MINT(84)+4 
        IF(MSTP(71).GE.1.AND.ISUB.NE.95.AND.K(IPU3,1).GT.0.AND. 
     &  K(IPU3,1).LE.10.AND.K(IPU4,1).GT.0.AND.K(IPU4,1).LE.10) THEN    
          QMAX=SQRT(PARP(71)*VINT(52))  
          IF(ISUB.EQ.5) QMAX=SQRT(PMAS(23,1)**2)    
          IF(ISUB.EQ.8) QMAX=SQRT(PMAS(24,1)**2)    
          CALL LUSHOW_HIJING(IPU3,IPU4,QMAX)   
        ENDIF   
    
C...Sum up transverse and longitudinal momenta. 
        IF(IOVL.EQ.1) THEN  
          PARI(65)=2.*PARI(17)  
          DO 130 I=MSTP(126)+1,N    
          IF(K(I,1).LE.0.OR.K(I,1).GT.10) GOTO 130  
          PT=SQRT(P(I,1)**2+P(I,2)**2)  
          PARI(69)=PARI(69)+PT  
          IF(I.LE.NSAV1.OR.I.GT.NSAV3) PARI(66)=PARI(66)+PT 
          IF(I.GT.NSAV1.AND.I.LE.NSAV2) PARI(68)=PARI(68)+PT    
  130     CONTINUE  
          PARI(67)=PARI(68) 
          PARI(71)=VINT(151)    
          PARI(72)=VINT(152)    
          PARI(73)=VINT(151)    
          PARI(74)=VINT(152)    
        ENDIF   
    
C...Decay of final state resonances.    
        IF(MSTP(41).GE.1.AND.ISUB.NE.95) CALL PYRESD_HIJING    
    
      ELSE  
C...Diffractive and elastic scattering. 
        CALL PYDIFF_HIJING 
        IF(IOVL.EQ.1) THEN  
          PARI(65)=2.*PARI(17)  
          PARI(66)=PARI(65) 
          PARI(69)=PARI(65) 
        ENDIF   
      ENDIF 
    
C...Recalculate energies from momenta and masses (if desired).  
      IF(MSTP(113).GE.1) THEN   
        DO 140 I=MINT(83)+1,N   
  140   IF(K(I,1).GT.0.AND.K(I,1).LE.10) P(I,4)=SQRT(P(I,1)**2+ 
     &  P(I,2)**2+P(I,3)**2+P(I,5)**2)  
      ENDIF 
    
C...Rearrange partons along strings, check invariant mass cuts. 
      MSTU(28)=0    
      CALL LUPREP_HIJING(MINT(84)+1)   
      IF(MSTP(112).EQ.1.AND.MSTU(28).EQ.3) GOTO 100 
      IF(MSTP(125).EQ.0.OR.MSTP(125).EQ.1) THEN 
        DO 150 I=MINT(84)+1,N   
        IF(K(I,2).NE.94) GOTO 150   
        K(I+1,3)=MOD(K(I+1,4)/MSTU(5),MSTU(5))  
        K(I+2,3)=MOD(K(I+2,4)/MSTU(5),MSTU(5))  
  150   CONTINUE    
        CALL LUEDIT_HIJING(12) 
        CALL LUEDIT_HIJING(14) 
        IF(MSTP(125).EQ.0) CALL LUEDIT_HIJING(15)  
        IF(MSTP(125).EQ.0) MINT(4)=0    
      ENDIF 
    
C...Introduce separators between sections in LULIST_HIJING event listing.  
      IF(IOVL.EQ.1.AND.MSTP(125).LE.0) THEN 
        MSTU(70)=1  
        MSTU(71)=N  
      ELSEIF(IOVL.EQ.1) THEN    
        MSTU(70)=3  
        MSTU(71)=2  
        MSTU(72)=MINT(4)    
        MSTU(73)=N  
      ENDIF 
    
C...Perform hadronization (if desired). 
      IF(MSTP(111).GE.1) CALL LUEXEC_HIJING    
      IF(MSTP(125).EQ.0.OR.MSTP(125).EQ.1) CALL LUEDIT_HIJING(14)  
    
C...Calculate Monte Carlo estimates of cross-sections.  
  160 IF(IOVL.EQ.1) THEN    
        IF(MSTP(111).NE.-1) NGEN(ISUB,3)=NGEN(ISUB,3)+1 
        NGEN(0,3)=NGEN(0,3)+1   
        XSEC(0,3)=0.    
        DO 170 I=1,200  
        IF(I.EQ.96) THEN    
          XSEC(I,3)=0.  
        ELSEIF(MSUB(95).EQ.1.AND.(I.EQ.11.OR.I.EQ.12.OR.I.EQ.13.OR. 
     &  I.EQ.28.OR.I.EQ.53.OR.I.EQ.68)) THEN    
          XSEC(I,3)=XSEC(96,2)*NGEN(I,3)/MAX(1.,FLOAT(NGEN(96,1))*  
     &    FLOAT(NGEN(96,2)))    
        ELSEIF(NGEN(I,1).EQ.0) THEN 
          XSEC(I,3)=0.  
        ELSEIF(NGEN(I,2).EQ.0) THEN 
          XSEC(I,3)=XSEC(I,2)*NGEN(0,3)/(FLOAT(NGEN(I,1))*  
     &    FLOAT(NGEN(0,2))) 
        ELSE    
          XSEC(I,3)=XSEC(I,2)*NGEN(I,3)/(FLOAT(NGEN(I,1))*  
     &    FLOAT(NGEN(I,2))) 
        ENDIF   
  170   XSEC(0,3)=XSEC(0,3)+XSEC(I,3)   
        IF(MSUB(95).EQ.1) THEN  
          NGENS=NGEN(91,3)+NGEN(92,3)+NGEN(93,3)+NGEN(94,3)+NGEN(95,3)  
          XSECS=XSEC(91,3)+XSEC(92,3)+XSEC(93,3)+XSEC(94,3)+XSEC(95,3)  
          XMAXS=XSEC(95,1)  
          IF(MSUB(91).EQ.1) XMAXS=XMAXS+XSEC(91,1)  
          IF(MSUB(92).EQ.1) XMAXS=XMAXS+XSEC(92,1)  
          IF(MSUB(93).EQ.1) XMAXS=XMAXS+XSEC(93,1)  
          IF(MSUB(94).EQ.1) XMAXS=XMAXS+XSEC(94,1)  
          FAC=1.    
          IF(NGENS.LT.NGEN(0,3)) FAC=(XMAXS-XSECS)/(XSEC(0,3)-XSECS)    
          XSEC(11,3)=FAC*XSEC(11,3) 
          XSEC(12,3)=FAC*XSEC(12,3) 
          XSEC(13,3)=FAC*XSEC(13,3) 
          XSEC(28,3)=FAC*XSEC(28,3) 
          XSEC(53,3)=FAC*XSEC(53,3) 
          XSEC(68,3)=FAC*XSEC(68,3) 
          XSEC(0,3)=XSEC(91,3)+XSEC(92,3)+XSEC(93,3)+XSEC(94,3)+    
     &    XSEC(95,1)    
        ENDIF   
    
C...Store final information.    
        MINT(5)=MINT(5)+1   
        MSTI(3)=MINT(3) 
        MSTI(4)=MINT(4) 
        MSTI(5)=MINT(5) 
        MSTI(6)=MINT(6) 
        MSTI(7)=MINT(7) 
        MSTI(8)=MINT(8) 
        MSTI(13)=MINT(13)   
        MSTI(14)=MINT(14)   
        MSTI(21)=MINT(21)   
        MSTI(22)=MINT(22)   
        MSTI(23)=MINT(23)   
        MSTI(24)=MINT(24)   
        MSTI(25)=MINT(25)   
        MSTI(26)=MINT(26)   
        MSTI(31)=MINT(31)   
        PARI(1)=XSEC(0,3)   
        PARI(2)=XSEC(0,3)/MINT(5)   
        PARI(31)=VINT(141)  
        PARI(32)=VINT(142)  
        IF(ISUB.NE.95.AND.MINT(7)*MINT(8).NE.0) THEN    
          PARI(42)=2.*VINT(47)/VINT(1)  
          DO 180 IS=7,8 
          PARI(36+IS)=P(MINT(IS),3)/VINT(1) 
          PARI(38+IS)=P(MINT(IS),4)/VINT(1) 
          I=MINT(IS)    
          PR=MAX(1E-20,P(I,5)**2+P(I,1)**2+P(I,2)**2)   
          PARI(40+IS)=SIGN(LOG(MIN((SQRT(PR+P(I,3)**2)+ABS(P(I,3)))/    
     &    SQRT(PR),1E20)),P(I,3))   
          PR=MAX(1E-20,P(I,1)**2+P(I,2)**2) 
          PARI(42+IS)=SIGN(LOG(MIN((SQRT(PR+P(I,3)**2)+ABS(P(I,3)))/    
     &    SQRT(PR),1E20)),P(I,3))   
          PARI(44+IS)=P(I,3)/SQRT(P(I,1)**2+P(I,2)**2+P(I,3)**2)    
          PARI(46+IS)=ULANGL_HIJING(P(I,3),SQRT(P(I,1)**2+P(I,2)**2))  
          PARI(48+IS)=ULANGL_HIJING(P(I,1),P(I,2)) 
  180     CONTINUE  
        ENDIF   
        PARI(61)=VINT(148)  
        IF(ISET(ISUB).EQ.1.OR.ISET(ISUB).EQ.3) THEN 
          MSTU(161)=MINT(21)    
          MSTU(162)=0   
        ELSE    
          MSTU(161)=MINT(21)    
          MSTU(162)=MINT(22)    
        ENDIF   
      ENDIF 
    
C...Prepare to go to next overlayed event.  
      MSTI(41)=IOVL 
      IF(IOVL.GE.2.AND.IOVL.LE.10) MSTI(40+IOVL)=ISUB   
      IF(MSTU(70).LT.10) THEN   
        MSTU(70)=MSTU(70)+1 
        MSTU(70+MSTU(70))=N 
      ENDIF 
      MINT(83)=N    
      MINT(84)=N+MSTP(126)  
  190 CONTINUE  
    
C...Information on overlayed events.    
      IF(MSTP(131).EQ.1.AND.MSTP(133).GE.1) THEN    
        PARI(91)=VINT(132)  
        PARI(92)=VINT(133)  
        PARI(93)=VINT(134)  
        IF(MSTP(133).EQ.2) PARI(93)=PARI(93)*XSEC(0,3)/VINT(131)    
      ENDIF 
    
C...Transform to the desired coordinate frame.  
  200 CALL PYFRAM_HIJING(MSTP(124))    
    
      RETURN    
      END   
