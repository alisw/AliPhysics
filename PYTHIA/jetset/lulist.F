 
C********************************************************************* 
 
      SUBROUTINE LULIST(MLIST) 
 
C...Purpose: to give program heading, or list an event, or particle 
C...data, or current parameter values. 
      COMMON/LUJETS/N,K(4000,5),P(4000,5),V(4000,5) 
      COMMON/LUDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200) 
      COMMON/LUDAT2/KCHG(500,3),PMAS(500,4),PARF(2000),VCKM(4,4) 
      COMMON/LUDAT3/MDCY(500,3),MDME(2000,2),BRAT(2000),KFDP(2000,5) 
      SAVE /LUJETS/,/LUDAT1/,/LUDAT2/,/LUDAT3/ 
      CHARACTER CHAP*16,CHAC*16,CHAN*16,CHAD(5)*16,CHDL(7)*4 
      DIMENSION PS(6) 
      DATA CHDL/'(())',' ','()','!!','<>','==','(==)'/ 
 
C...Initialization printout: version number and date of last change. 
      IF(MLIST.EQ.0.OR.MSTU(12).EQ.1) THEN 
        CALL LULOGO 
        MSTU(12)=0 
        IF(MLIST.EQ.0) RETURN 
      ENDIF 
 
C...List event data, including additional lines after N. 
      IF(MLIST.GE.1.AND.MLIST.LE.3) THEN 
        IF(MLIST.EQ.1) WRITE(MSTU(11),5100) 
        IF(MLIST.EQ.2) WRITE(MSTU(11),5200) 
        IF(MLIST.EQ.3) WRITE(MSTU(11),5300) 
        LMX=12 
        IF(MLIST.GE.2) LMX=16 
        ISTR=0 
        IMAX=N 
        IF(MSTU(2).GT.0) IMAX=MSTU(2) 
        DO 120 I=MAX(1,MSTU(1)),MAX(IMAX,N+MAX(0,MSTU(3))) 
        IF((I.GT.IMAX.AND.I.LE.N).OR.K(I,1).LT.0) GOTO 120 
 
C...Get particle name, pad it and check it is not too long. 
        CALL LUNAME(K(I,2),CHAP) 
        LEN=0 
        DO 100 LEM=1,16 
        IF(CHAP(LEM:LEM).NE.' ') LEN=LEM 
  100   CONTINUE 
        MDL=(K(I,1)+19)/10 
        LDL=0 
        IF(MDL.EQ.2.OR.MDL.GE.8) THEN 
          CHAC=CHAP 
          IF(LEN.GT.LMX) CHAC(LMX:LMX)='?' 
        ELSE 
          LDL=1 
          IF(MDL.EQ.1.OR.MDL.EQ.7) LDL=2 
          IF(LEN.EQ.0) THEN 
            CHAC=CHDL(MDL)(1:2*LDL)//' ' 
          ELSE 
            CHAC=CHDL(MDL)(1:LDL)//CHAP(1:MIN(LEN,LMX-2*LDL))// 
     &      CHDL(MDL)(LDL+1:2*LDL)//' ' 
            IF(LEN+2*LDL.GT.LMX) CHAC(LMX:LMX)='?' 
          ENDIF 
        ENDIF 
 
C...Add information on string connection. 
        IF(K(I,1).EQ.1.OR.K(I,1).EQ.2.OR.K(I,1).EQ.11.OR.K(I,1).EQ.12) 
     &  THEN 
          KC=LUCOMP(K(I,2)) 
          KCC=0 
          IF(KC.NE.0) KCC=KCHG(KC,2) 
          IF(IABS(K(I,2)).EQ.39) THEN 
            IF(LEN+2*LDL+3.LE.LMX) CHAC(LMX-1:LMX-1)='X' 
          ELSEIF(KCC.NE.0.AND.ISTR.EQ.0) THEN 
            ISTR=1 
            IF(LEN+2*LDL+3.LE.LMX) CHAC(LMX-1:LMX-1)='A' 
          ELSEIF(KCC.NE.0.AND.(K(I,1).EQ.2.OR.K(I,1).EQ.12)) THEN 
            IF(LEN+2*LDL+3.LE.LMX) CHAC(LMX-1:LMX-1)='I' 
          ELSEIF(KCC.NE.0) THEN 
            ISTR=0 
            IF(LEN+2*LDL+3.LE.LMX) CHAC(LMX-1:LMX-1)='V' 
          ENDIF 
        ENDIF 
 
C...Write data for particle/jet. 
        IF(MLIST.EQ.1.AND.ABS(P(I,4)).LT.9999.) THEN 
          WRITE(MSTU(11),5400) I,CHAC(1:12),(K(I,J1),J1=1,3), 
     &    (P(I,J2),J2=1,5) 
        ELSEIF(MLIST.EQ.1.AND.ABS(P(I,4)).LT.99999.) THEN 
          WRITE(MSTU(11),5500) I,CHAC(1:12),(K(I,J1),J1=1,3), 
     &    (P(I,J2),J2=1,5) 
        ELSEIF(MLIST.EQ.1) THEN 
          WRITE(MSTU(11),5600) I,CHAC(1:12),(K(I,J1),J1=1,3), 
     &    (P(I,J2),J2=1,5) 
        ELSEIF(MSTU(5).EQ.10000.AND.(K(I,1).EQ.3.OR.K(I,1).EQ.13.OR. 
     &  K(I,1).EQ.14)) THEN 
          WRITE(MSTU(11),5700) I,CHAC,(K(I,J1),J1=1,3), 
     &    K(I,4)/100000000,MOD(K(I,4)/10000,10000),MOD(K(I,4),10000), 
     &    K(I,5)/100000000,MOD(K(I,5)/10000,10000),MOD(K(I,5),10000), 
     &    (P(I,J2),J2=1,5) 
        ELSE 
          WRITE(MSTU(11),5800) I,CHAC,(K(I,J1),J1=1,5),(P(I,J2),J2=1,5) 
        ENDIF 
        IF(MLIST.EQ.3) WRITE(MSTU(11),5900) (V(I,J),J=1,5) 
 
C...Insert extra separator lines specified by user. 
        IF(MSTU(70).GE.1) THEN 
          ISEP=0 
          DO 110 J=1,MIN(10,MSTU(70)) 
          IF(I.EQ.MSTU(70+J)) ISEP=1 
  110     CONTINUE 
          IF(ISEP.EQ.1.AND.MLIST.EQ.1) WRITE(MSTU(11),6000) 
          IF(ISEP.EQ.1.AND.MLIST.GE.2) WRITE(MSTU(11),6100) 
        ENDIF 
  120   CONTINUE 
 
C...Sum of charges and momenta. 
        DO 130 J=1,6 
        PS(J)=PLU(0,J) 
  130   CONTINUE 
        IF(MLIST.EQ.1.AND.ABS(PS(4)).LT.9999.) THEN 
          WRITE(MSTU(11),6200) PS(6),(PS(J),J=1,5) 
        ELSEIF(MLIST.EQ.1.AND.ABS(PS(4)).LT.99999.) THEN 
          WRITE(MSTU(11),6300) PS(6),(PS(J),J=1,5) 
        ELSEIF(MLIST.EQ.1) THEN 
          WRITE(MSTU(11),6400) PS(6),(PS(J),J=1,5) 
        ELSE 
          WRITE(MSTU(11),6500) PS(6),(PS(J),J=1,5) 
        ENDIF 
 
C...Give simple list of KF codes defined in program. 
      ELSEIF(MLIST.EQ.11) THEN 
        WRITE(MSTU(11),6600) 
        DO 140 KF=1,40 
        CALL LUNAME(KF,CHAP) 
        CALL LUNAME(-KF,CHAN) 
        IF(CHAP.NE.' '.AND.CHAN.EQ.' ') WRITE(MSTU(11),6700) KF,CHAP 
        IF(CHAN.NE.' ') WRITE(MSTU(11),6700) KF,CHAP,-KF,CHAN 
  140   CONTINUE 
        DO 170 KFLS=1,3,2 
        DO 160 KFLA=1,8 
        DO 150 KFLB=1,KFLA-(3-KFLS)/2 
        KF=1000*KFLA+100*KFLB+KFLS 
        CALL LUNAME(KF,CHAP) 
        CALL LUNAME(-KF,CHAN) 
        WRITE(MSTU(11),6700) KF,CHAP,-KF,CHAN 
  150   CONTINUE 
  160   CONTINUE 
  170   CONTINUE 
        KF=130 
        CALL LUNAME(KF,CHAP) 
        WRITE(MSTU(11),6700) KF,CHAP 
        KF=310 
        CALL LUNAME(KF,CHAP) 
        WRITE(MSTU(11),6700) KF,CHAP 
        DO 200 KMUL=0,5 
        KFLS=3 
        IF(KMUL.EQ.0.OR.KMUL.EQ.3) KFLS=1 
        IF(KMUL.EQ.5) KFLS=5 
        KFLR=0 
        IF(KMUL.EQ.2.OR.KMUL.EQ.3) KFLR=1 
        IF(KMUL.EQ.4) KFLR=2 
        DO 190 KFLB=1,8 
        DO 180 KFLC=1,KFLB-1 
        KF=10000*KFLR+100*KFLB+10*KFLC+KFLS 
        CALL LUNAME(KF,CHAP) 
        CALL LUNAME(-KF,CHAN) 
        WRITE(MSTU(11),6700) KF,CHAP,-KF,CHAN 
  180   CONTINUE 
        KF=10000*KFLR+110*KFLB+KFLS 
        CALL LUNAME(KF,CHAP) 
        WRITE(MSTU(11),6700) KF,CHAP 
  190   CONTINUE 
  200 CONTINUE 
        KF=30443 
        CALL LUNAME(KF,CHAP) 
        WRITE(MSTU(11),6700) KF,CHAP 
        KF=30553 
        CALL LUNAME(KF,CHAP) 
        WRITE(MSTU(11),6700) KF,CHAP 
        DO 240 KFLSP=1,3 
        KFLS=2+2*(KFLSP/3) 
        DO 230 KFLA=1,8 
        DO 220 KFLB=1,KFLA 
        DO 210 KFLC=1,KFLB 
        IF(KFLSP.EQ.1.AND.(KFLA.EQ.KFLB.OR.KFLB.EQ.KFLC)) GOTO 210 
        IF(KFLSP.EQ.2.AND.KFLA.EQ.KFLC) GOTO 210 
        IF(KFLSP.EQ.1) KF=1000*KFLA+100*KFLC+10*KFLB+KFLS 
        IF(KFLSP.GE.2) KF=1000*KFLA+100*KFLB+10*KFLC+KFLS 
        CALL LUNAME(KF,CHAP) 
        CALL LUNAME(-KF,CHAN) 
        WRITE(MSTU(11),6700) KF,CHAP,-KF,CHAN 
  210   CONTINUE 
  220   CONTINUE 
  230   CONTINUE 
  240   CONTINUE 
 
C...List parton/particle data table. Check whether to be listed. 
      ELSEIF(MLIST.EQ.12) THEN 
        WRITE(MSTU(11),6800) 
        MSTJ24=MSTJ(24) 
        MSTJ(24)=0 
        KFMAX=30553 
        IF(MSTU(2).NE.0) KFMAX=MSTU(2) 
        DO 270 KF=MAX(1,MSTU(1)),KFMAX 
        KC=LUCOMP(KF) 
        IF(KC.EQ.0) GOTO 270 
        IF(MSTU(14).EQ.0.AND.KF.GT.100.AND.KC.LE.100) GOTO 270 
        IF(MSTU(14).GT.0.AND.KF.GT.100.AND.MAX(MOD(KF/1000,10), 
     &  MOD(KF/100,10)).GT.MSTU(14)) GOTO 270 
        IF(MSTU(14).GT.0.AND.KF.GT.100.AND.KC.EQ.90) GOTO 270 
 
C...Find particle name and mass. Print information. 
        CALL LUNAME(KF,CHAP) 
        IF(KF.LE.100.AND.CHAP.EQ.' '.AND.MDCY(KC,2).EQ.0) GOTO 270 
        CALL LUNAME(-KF,CHAN) 
        PM=ULMASS(KF) 
        WRITE(MSTU(11),6900) KF,KC,CHAP,CHAN,KCHG(KC,1),KCHG(KC,2), 
     &  KCHG(KC,3),PM,PMAS(KC,2),PMAS(KC,3),PMAS(KC,4),MDCY(KC,1) 
 
C...Particle decay: channel number, branching ration, matrix element, 
C...decay products. 
        IF(KF.GT.100.AND.KC.LE.100) GOTO 270 
        DO 260 IDC=MDCY(KC,2),MDCY(KC,2)+MDCY(KC,3)-1 
        DO 250 J=1,5 
        CALL LUNAME(KFDP(IDC,J),CHAD(J)) 
  250   CONTINUE 
        WRITE(MSTU(11),7000) IDC,MDME(IDC,1),MDME(IDC,2),BRAT(IDC), 
     &  (CHAD(J),J=1,5) 
  260   CONTINUE 
  270   CONTINUE 
        MSTJ(24)=MSTJ24 
 
C...List parameter value table. 
      ELSEIF(MLIST.EQ.13) THEN 
        WRITE(MSTU(11),7100) 
        DO 280 I=1,200 
        WRITE(MSTU(11),7200) I,MSTU(I),PARU(I),MSTJ(I),PARJ(I),PARF(I) 
  280   CONTINUE 
      ENDIF 
 
C...Format statements for output on unit MSTU(11) (by default 6). 
 5100 FORMAT(///28X,'Event listing (summary)'//4X,'I  particle/jet KS', 
     &5X,'KF orig    p_x      p_y      p_z       E        m'/) 
 5200 FORMAT(///28X,'Event listing (standard)'//4X,'I  particle/jet', 
     &'  K(I,1)   K(I,2) K(I,3)     K(I,4)      K(I,5)       P(I,1)', 
     &'       P(I,2)       P(I,3)       P(I,4)       P(I,5)'/) 
 5300 FORMAT(///28X,'Event listing (with vertices)'//4X,'I  particle/j', 
     &'et  K(I,1)   K(I,2) K(I,3)     K(I,4)      K(I,5)       P(I,1)', 
     &'       P(I,2)       P(I,3)       P(I,4)       P(I,5)'/73X, 
     &'V(I,1)       V(I,2)       V(I,3)       V(I,4)       V(I,5)'/) 
 5400 FORMAT(1X,I4,2X,A12,1X,I2,1X,I6,1X,I4,5F9.3) 
 5500 FORMAT(1X,I4,2X,A12,1X,I2,1X,I6,1X,I4,5F9.2) 
 5600 FORMAT(1X,I4,2X,A12,1X,I2,1X,I6,1X,I4,5F9.1) 
 5700 FORMAT(1X,I4,2X,A16,1X,I3,1X,I8,2X,I4,2(3X,I1,2I4),5F13.5) 
 5800 FORMAT(1X,I4,2X,A16,1X,I3,1X,I8,2X,I4,2(3X,I9),5F13.5) 
 5900 FORMAT(66X,5(1X,F12.3)) 
 6000 FORMAT(1X,78('=')) 
 6100 FORMAT(1X,130('=')) 
 6200 FORMAT(19X,'sum:',F6.2,5X,5F9.3) 
 6300 FORMAT(19X,'sum:',F6.2,5X,5F9.2) 
 6400 FORMAT(19X,'sum:',F6.2,5X,5F9.1) 
 6500 FORMAT(19X,'sum charge:',F6.2,3X,'sum momentum and inv. mass:', 
     &5F13.5) 
 6600 FORMAT(///20X,'List of KF codes in program'/) 
 6700 FORMAT(4X,I6,4X,A16,6X,I6,4X,A16) 
 6800 FORMAT(///30X,'Particle/parton data table'//5X,'KF',5X,'KC',4X, 
     &'particle',8X,'antiparticle',6X,'chg  col  anti',8X,'mass',7X, 
     &'width',7X,'w-cut',5X,'lifetime',1X,'decay'/11X,'IDC',1X,'on/off', 
     &1X,'ME',3X,'Br.rat.',4X,'decay products') 
 6900 FORMAT(/1X,I6,3X,I4,4X,A16,A16,3I5,1X,F12.5,2(1X,F11.5), 
     &2X,F12.5,3X,I2) 
 7000 FORMAT(10X,I4,2X,I3,2X,I3,2X,F8.5,4X,5A16) 
 7100 FORMAT(///20X,'Parameter value table'//4X,'I',3X,'MSTU(I)', 
     &8X,'PARU(I)',3X,'MSTJ(I)',8X,'PARJ(I)',8X,'PARF(I)') 
 7200 FORMAT(1X,I4,1X,I9,1X,F14.5,1X,I9,1X,F14.5,1X,F14.5) 
 
      RETURN 
      END 
