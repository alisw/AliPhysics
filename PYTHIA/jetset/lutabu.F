 
C********************************************************************* 
 
      SUBROUTINE LUTABU(MTABU) 
 
C...Purpose: to evaluate various properties of an event, with 
C...statistics accumulated during the course of the run and 
C...printed at the end. 
      COMMON/LUJETS/N,K(4000,5),P(4000,5),V(4000,5) 
      COMMON/LUDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200) 
      COMMON/LUDAT2/KCHG(500,3),PMAS(500,4),PARF(2000),VCKM(4,4) 
      COMMON/LUDAT3/MDCY(500,3),MDME(2000,2),BRAT(2000),KFDP(2000,5) 
      SAVE /LUJETS/,/LUDAT1/,/LUDAT2/,/LUDAT3/ 
      DIMENSION KFIS(100,2),NPIS(100,0:10),KFFS(400),NPFS(400,4), 
     &FEVFM(10,4),FM1FM(3,10,4),FM2FM(3,10,4),FMOMA(4),FMOMS(4), 
     &FEVEE(50),FE1EC(50),FE2EC(50),FE1EA(25),FE2EA(25), 
     &KFDM(8),KFDC(200,0:8),NPDC(200) 
      SAVE NEVIS,NKFIS,KFIS,NPIS,NEVFS,NPRFS,NFIFS,NCHFS,NKFFS, 
     &KFFS,NPFS,NEVFM,NMUFM,FM1FM,FM2FM,NEVEE,FE1EC,FE2EC,FE1EA, 
     &FE2EA,NEVDC,NKFDC,NREDC,KFDC,NPDC 
      CHARACTER CHAU*16,CHIS(2)*12,CHDC(8)*12 
      DATA NEVIS/0/,NKFIS/0/,NEVFS/0/,NPRFS/0/,NFIFS/0/,NCHFS/0/, 
     &NKFFS/0/,NEVFM/0/,NMUFM/0/,FM1FM/120*0./,FM2FM/120*0./, 
     &NEVEE/0/,FE1EC/50*0./,FE2EC/50*0./,FE1EA/25*0./,FE2EA/25*0./, 
     &NEVDC/0/,NKFDC/0/,NREDC/0/ 
 
C...Reset statistics on initial parton state. 
      IF(MTABU.EQ.10) THEN 
        NEVIS=0 
        NKFIS=0 
 
C...Identify and order flavour content of initial state. 
      ELSEIF(MTABU.EQ.11) THEN 
        NEVIS=NEVIS+1 
        KFM1=2*IABS(MSTU(161)) 
        IF(MSTU(161).GT.0) KFM1=KFM1-1 
        KFM2=2*IABS(MSTU(162)) 
        IF(MSTU(162).GT.0) KFM2=KFM2-1 
        KFMN=MIN(KFM1,KFM2) 
        KFMX=MAX(KFM1,KFM2) 
        DO 100 I=1,NKFIS 
        IF(KFMN.EQ.KFIS(I,1).AND.KFMX.EQ.KFIS(I,2)) THEN 
          IKFIS=-I 
          GOTO 110 
        ELSEIF(KFMN.LT.KFIS(I,1).OR.(KFMN.EQ.KFIS(I,1).AND. 
     &  KFMX.LT.KFIS(I,2))) THEN 
          IKFIS=I 
          GOTO 110 
        ENDIF 
  100   CONTINUE 
        IKFIS=NKFIS+1 
  110   IF(IKFIS.LT.0) THEN 
          IKFIS=-IKFIS 
        ELSE 
          IF(NKFIS.GE.100) RETURN 
          DO 130 I=NKFIS,IKFIS,-1 
          KFIS(I+1,1)=KFIS(I,1) 
          KFIS(I+1,2)=KFIS(I,2) 
          DO 120 J=0,10 
          NPIS(I+1,J)=NPIS(I,J) 
  120     CONTINUE 
  130   CONTINUE 
          NKFIS=NKFIS+1 
          KFIS(IKFIS,1)=KFMN 
          KFIS(IKFIS,2)=KFMX 
          DO 140 J=0,10 
          NPIS(IKFIS,J)=0 
  140     CONTINUE 
        ENDIF 
        NPIS(IKFIS,0)=NPIS(IKFIS,0)+1 
 
C...Count number of partons in initial state. 
        NP=0 
        DO 160 I=1,N 
        IF(K(I,1).LE.0.OR.K(I,1).GT.12) THEN 
        ELSEIF(IABS(K(I,2)).GT.80.AND.IABS(K(I,2)).LE.100) THEN 
        ELSEIF(IABS(K(I,2)).GT.100.AND.MOD(IABS(K(I,2))/10,10).NE.0) 
     &  THEN 
        ELSE 
          IM=I 
  150     IM=K(IM,3) 
          IF(IM.LE.0.OR.IM.GT.N) THEN 
            NP=NP+1 
          ELSEIF(K(IM,1).LE.0.OR.K(IM,1).GT.20) THEN 
            NP=NP+1 
          ELSEIF(IABS(K(IM,2)).GT.80.AND.IABS(K(IM,2)).LE.100) THEN 
          ELSEIF(IABS(K(IM,2)).GT.100.AND.MOD(IABS(K(IM,2))/10,10).NE.0) 
     &    THEN 
          ELSE 
            GOTO 150 
          ENDIF 
        ENDIF 
  160   CONTINUE 
        NPCO=MAX(NP,1) 
        IF(NP.GE.6) NPCO=6 
        IF(NP.GE.8) NPCO=7 
        IF(NP.GE.11) NPCO=8 
        IF(NP.GE.16) NPCO=9 
        IF(NP.GE.26) NPCO=10 
        NPIS(IKFIS,NPCO)=NPIS(IKFIS,NPCO)+1 
        MSTU(62)=NP 
 
C...Write statistics on initial parton state. 
      ELSEIF(MTABU.EQ.12) THEN 
        FAC=1./MAX(1,NEVIS) 
        WRITE(MSTU(11),5000) NEVIS 
        DO 170 I=1,NKFIS 
        KFMN=KFIS(I,1) 
        IF(KFMN.EQ.0) KFMN=KFIS(I,2) 
        KFM1=(KFMN+1)/2 
        IF(2*KFM1.EQ.KFMN) KFM1=-KFM1 
        CALL LUNAME(KFM1,CHAU) 
        CHIS(1)=CHAU(1:12) 
        IF(CHAU(13:13).NE.' ') CHIS(1)(12:12)='?' 
        KFMX=KFIS(I,2) 
        IF(KFIS(I,1).EQ.0) KFMX=0 
        KFM2=(KFMX+1)/2 
        IF(2*KFM2.EQ.KFMX) KFM2=-KFM2 
        CALL LUNAME(KFM2,CHAU) 
        CHIS(2)=CHAU(1:12) 
        IF(CHAU(13:13).NE.' ') CHIS(2)(12:12)='?' 
        WRITE(MSTU(11),5100) CHIS(1),CHIS(2),FAC*NPIS(I,0), 
     &  (NPIS(I,J)/FLOAT(NPIS(I,0)),J=1,10) 
  170   CONTINUE 
 
C...Copy statistics on initial parton state into /LUJETS/. 
      ELSEIF(MTABU.EQ.13) THEN 
        FAC=1./MAX(1,NEVIS) 
        DO 190 I=1,NKFIS 
        KFMN=KFIS(I,1) 
        IF(KFMN.EQ.0) KFMN=KFIS(I,2) 
        KFM1=(KFMN+1)/2 
        IF(2*KFM1.EQ.KFMN) KFM1=-KFM1 
        KFMX=KFIS(I,2) 
        IF(KFIS(I,1).EQ.0) KFMX=0 
        KFM2=(KFMX+1)/2 
        IF(2*KFM2.EQ.KFMX) KFM2=-KFM2 
        K(I,1)=32 
        K(I,2)=99 
        K(I,3)=KFM1 
        K(I,4)=KFM2 
        K(I,5)=NPIS(I,0) 
        DO 180 J=1,5 
        P(I,J)=FAC*NPIS(I,J) 
        V(I,J)=FAC*NPIS(I,J+5) 
  180   CONTINUE 
  190   CONTINUE 
        N=NKFIS 
        DO 200 J=1,5 
        K(N+1,J)=0 
        P(N+1,J)=0. 
        V(N+1,J)=0. 
  200   CONTINUE 
        K(N+1,1)=32 
        K(N+1,2)=99 
        K(N+1,5)=NEVIS 
        MSTU(3)=1 
 
C...Reset statistics on number of particles/partons. 
      ELSEIF(MTABU.EQ.20) THEN 
        NEVFS=0 
        NPRFS=0 
        NFIFS=0 
        NCHFS=0 
        NKFFS=0 
 
C...Identify whether particle/parton is primary or not. 
      ELSEIF(MTABU.EQ.21) THEN 
        NEVFS=NEVFS+1 
        MSTU(62)=0 
        DO 260 I=1,N 
        IF(K(I,1).LE.0.OR.K(I,1).GT.20.OR.K(I,1).EQ.13) GOTO 260 
        MSTU(62)=MSTU(62)+1 
        KC=LUCOMP(K(I,2)) 
        MPRI=0 
        IF(K(I,3).LE.0.OR.K(I,3).GT.N) THEN 
          MPRI=1 
        ELSEIF(K(K(I,3),1).LE.0.OR.K(K(I,3),1).GT.20) THEN 
          MPRI=1 
        ELSEIF(K(K(I,3),2).GE.91.AND.K(K(I,3),2).LE.93) THEN 
          MPRI=1 
        ELSEIF(KC.EQ.0) THEN 
        ELSEIF(K(K(I,3),1).EQ.13) THEN 
          IM=K(K(I,3),3) 
          IF(IM.LE.0.OR.IM.GT.N) THEN 
            MPRI=1 
          ELSEIF(K(IM,1).LE.0.OR.K(IM,1).GT.20) THEN 
            MPRI=1 
          ENDIF 
        ELSEIF(KCHG(KC,2).EQ.0) THEN 
          KCM=LUCOMP(K(K(I,3),2)) 
          IF(KCM.NE.0) THEN 
            IF(KCHG(KCM,2).NE.0) MPRI=1 
          ENDIF 
        ENDIF 
        IF(KC.NE.0.AND.MPRI.EQ.1) THEN 
          IF(KCHG(KC,2).EQ.0) NPRFS=NPRFS+1 
        ENDIF 
        IF(K(I,1).LE.10) THEN 
          NFIFS=NFIFS+1 
          IF(LUCHGE(K(I,2)).NE.0) NCHFS=NCHFS+1 
        ENDIF 
 
C...Fill statistics on number of particles/partons in event. 
        KFA=IABS(K(I,2)) 
        KFS=3-ISIGN(1,K(I,2))-MPRI 
        DO 210 IP=1,NKFFS 
        IF(KFA.EQ.KFFS(IP)) THEN 
          IKFFS=-IP 
          GOTO 220 
        ELSEIF(KFA.LT.KFFS(IP)) THEN 
          IKFFS=IP 
          GOTO 220 
        ENDIF 
  210   CONTINUE 
        IKFFS=NKFFS+1 
  220   IF(IKFFS.LT.0) THEN 
          IKFFS=-IKFFS 
        ELSE 
          IF(NKFFS.GE.400) RETURN 
          DO 240 IP=NKFFS,IKFFS,-1 
          KFFS(IP+1)=KFFS(IP) 
          DO 230 J=1,4 
          NPFS(IP+1,J)=NPFS(IP,J) 
  230     CONTINUE 
  240   CONTINUE 
          NKFFS=NKFFS+1 
          KFFS(IKFFS)=KFA 
          DO 250 J=1,4 
          NPFS(IKFFS,J)=0 
  250     CONTINUE 
        ENDIF 
        NPFS(IKFFS,KFS)=NPFS(IKFFS,KFS)+1 
  260   CONTINUE 
 
C...Write statistics on particle/parton composition of events. 
      ELSEIF(MTABU.EQ.22) THEN 
        FAC=1./MAX(1,NEVFS) 
        WRITE(MSTU(11),5200) NEVFS,FAC*NPRFS,FAC*NFIFS,FAC*NCHFS 
        DO 270 I=1,NKFFS 
        CALL LUNAME(KFFS(I),CHAU) 
        KC=LUCOMP(KFFS(I)) 
        MDCYF=0 
        IF(KC.NE.0) MDCYF=MDCY(KC,1) 
        WRITE(MSTU(11),5300) KFFS(I),CHAU,MDCYF,(FAC*NPFS(I,J),J=1,4), 
     &  FAC*(NPFS(I,1)+NPFS(I,2)+NPFS(I,3)+NPFS(I,4)) 
  270   CONTINUE 
 
C...Copy particle/parton composition information into /LUJETS/. 
      ELSEIF(MTABU.EQ.23) THEN 
        FAC=1./MAX(1,NEVFS) 
        DO 290 I=1,NKFFS 
        K(I,1)=32 
        K(I,2)=99 
        K(I,3)=KFFS(I) 
        K(I,4)=0 
        K(I,5)=NPFS(I,1)+NPFS(I,2)+NPFS(I,3)+NPFS(I,4) 
        DO 280 J=1,4 
        P(I,J)=FAC*NPFS(I,J) 
        V(I,J)=0. 
  280   CONTINUE 
        P(I,5)=FAC*K(I,5) 
        V(I,5)=0. 
  290   CONTINUE 
        N=NKFFS 
        DO 300 J=1,5 
        K(N+1,J)=0 
        P(N+1,J)=0. 
        V(N+1,J)=0. 
  300   CONTINUE 
        K(N+1,1)=32 
        K(N+1,2)=99 
        K(N+1,5)=NEVFS 
        P(N+1,1)=FAC*NPRFS 
        P(N+1,2)=FAC*NFIFS 
        P(N+1,3)=FAC*NCHFS 
        MSTU(3)=1 
 
C...Reset factorial moments statistics. 
      ELSEIF(MTABU.EQ.30) THEN 
        NEVFM=0 
        NMUFM=0 
        DO 330 IM=1,3 
        DO 320 IB=1,10 
        DO 310 IP=1,4 
        FM1FM(IM,IB,IP)=0. 
        FM2FM(IM,IB,IP)=0. 
  310   CONTINUE 
  320   CONTINUE 
  330   CONTINUE 
 
C...Find particles to include, with (pion,pseudo)rapidity and azimuth. 
      ELSEIF(MTABU.EQ.31) THEN 
        NEVFM=NEVFM+1 
        NLOW=N+MSTU(3) 
        NUPP=NLOW 
        DO 410 I=1,N 
        IF(K(I,1).LE.0.OR.K(I,1).GT.10) GOTO 410 
        IF(MSTU(41).GE.2) THEN 
          KC=LUCOMP(K(I,2)) 
          IF(KC.EQ.0.OR.KC.EQ.12.OR.KC.EQ.14.OR.KC.EQ.16.OR. 
     &    KC.EQ.18) GOTO 410 
          IF(MSTU(41).GE.3.AND.KCHG(KC,2).EQ.0.AND.LUCHGE(K(I,2)).EQ.0) 
     &    GOTO 410 
        ENDIF 
        PMR=0. 
        IF(MSTU(42).EQ.1.AND.K(I,2).NE.22) PMR=ULMASS(211) 
        IF(MSTU(42).GE.2) PMR=P(I,5) 
        PR=MAX(1E-20,PMR**2+P(I,1)**2+P(I,2)**2) 
        YETA=SIGN(LOG(MIN((SQRT(PR+P(I,3)**2)+ABS(P(I,3)))/SQRT(PR), 
     &  1E20)),P(I,3)) 
        IF(ABS(YETA).GT.PARU(57)) GOTO 410 
        PHI=ULANGL(P(I,1),P(I,2)) 
        IYETA=512.*(YETA+PARU(57))/(2.*PARU(57)) 
        IYETA=MAX(0,MIN(511,IYETA)) 
        IPHI=512.*(PHI+PARU(1))/PARU(2) 
        IPHI=MAX(0,MIN(511,IPHI)) 
        IYEP=0 
        DO 340 IB=0,9 
        IYEP=IYEP+4**IB*(2*MOD(IYETA/2**IB,2)+MOD(IPHI/2**IB,2)) 
  340   CONTINUE 
 
C...Order particles in (pseudo)rapidity and/or azimuth. 
        IF(NUPP.GT.MSTU(4)-5-MSTU(32)) THEN 
          CALL LUERRM(11,'(LUTABU:) no more memory left in LUJETS') 
          RETURN 
        ENDIF 
        NUPP=NUPP+1 
        IF(NUPP.EQ.NLOW+1) THEN 
          K(NUPP,1)=IYETA 
          K(NUPP,2)=IPHI 
          K(NUPP,3)=IYEP 
        ELSE 
          DO 350 I1=NUPP-1,NLOW+1,-1 
          IF(IYETA.GE.K(I1,1)) GOTO 360 
          K(I1+1,1)=K(I1,1) 
  350     CONTINUE 
  360     K(I1+1,1)=IYETA 
          DO 370 I1=NUPP-1,NLOW+1,-1 
          IF(IPHI.GE.K(I1,2)) GOTO 380 
          K(I1+1,2)=K(I1,2) 
  370     CONTINUE 
  380     K(I1+1,2)=IPHI 
          DO 390 I1=NUPP-1,NLOW+1,-1 
          IF(IYEP.GE.K(I1,3)) GOTO 400 
          K(I1+1,3)=K(I1,3) 
  390     CONTINUE 
  400     K(I1+1,3)=IYEP 
        ENDIF 
  410   CONTINUE 
        K(NUPP+1,1)=2**10 
        K(NUPP+1,2)=2**10 
        K(NUPP+1,3)=4**10 
 
C...Calculate sum of factorial moments in event. 
        DO 480 IM=1,3 
        DO 430 IB=1,10 
        DO 420 IP=1,4 
        FEVFM(IB,IP)=0. 
  420   CONTINUE 
  430   CONTINUE 
        DO 450 IB=1,10 
        IF(IM.LE.2) IBIN=2**(10-IB) 
        IF(IM.EQ.3) IBIN=4**(10-IB) 
        IAGR=K(NLOW+1,IM)/IBIN 
        NAGR=1 
        DO 440 I=NLOW+2,NUPP+1 
        ICUT=K(I,IM)/IBIN 
        IF(ICUT.EQ.IAGR) THEN 
          NAGR=NAGR+1 
        ELSE 
          IF(NAGR.EQ.1) THEN 
          ELSEIF(NAGR.EQ.2) THEN 
            FEVFM(IB,1)=FEVFM(IB,1)+2. 
          ELSEIF(NAGR.EQ.3) THEN 
            FEVFM(IB,1)=FEVFM(IB,1)+6. 
            FEVFM(IB,2)=FEVFM(IB,2)+6. 
          ELSEIF(NAGR.EQ.4) THEN 
            FEVFM(IB,1)=FEVFM(IB,1)+12. 
            FEVFM(IB,2)=FEVFM(IB,2)+24. 
            FEVFM(IB,3)=FEVFM(IB,3)+24. 
          ELSE 
            FEVFM(IB,1)=FEVFM(IB,1)+NAGR*(NAGR-1.) 
            FEVFM(IB,2)=FEVFM(IB,2)+NAGR*(NAGR-1.)*(NAGR-2.) 
            FEVFM(IB,3)=FEVFM(IB,3)+NAGR*(NAGR-1.)*(NAGR-2.)*(NAGR-3.) 
            FEVFM(IB,4)=FEVFM(IB,4)+NAGR*(NAGR-1.)*(NAGR-2.)*(NAGR-3.)* 
     &      (NAGR-4.) 
          ENDIF 
          IAGR=ICUT 
          NAGR=1 
        ENDIF 
  440   CONTINUE 
  450   CONTINUE 
 
C...Add results to total statistics. 
        DO 470 IB=10,1,-1 
        DO 460 IP=1,4 
        IF(FEVFM(1,IP).LT.0.5) THEN 
          FEVFM(IB,IP)=0. 
        ELSEIF(IM.LE.2) THEN 
          FEVFM(IB,IP)=2.**((IB-1)*IP)*FEVFM(IB,IP)/FEVFM(1,IP) 
        ELSE 
          FEVFM(IB,IP)=4.**((IB-1)*IP)*FEVFM(IB,IP)/FEVFM(1,IP) 
        ENDIF 
        FM1FM(IM,IB,IP)=FM1FM(IM,IB,IP)+FEVFM(IB,IP) 
        FM2FM(IM,IB,IP)=FM2FM(IM,IB,IP)+FEVFM(IB,IP)**2 
  460   CONTINUE 
  470   CONTINUE 
  480   CONTINUE 
        NMUFM=NMUFM+(NUPP-NLOW) 
        MSTU(62)=NUPP-NLOW 
 
C...Write accumulated statistics on factorial moments. 
      ELSEIF(MTABU.EQ.32) THEN 
        FAC=1./MAX(1,NEVFM) 
        IF(MSTU(42).LE.0) WRITE(MSTU(11),5400) NEVFM,'eta' 
        IF(MSTU(42).EQ.1) WRITE(MSTU(11),5400) NEVFM,'ypi' 
        IF(MSTU(42).GE.2) WRITE(MSTU(11),5400) NEVFM,'y  ' 
        DO 510 IM=1,3 
        WRITE(MSTU(11),5500) 
        DO 500 IB=1,10 
        BYETA=2.*PARU(57) 
        IF(IM.NE.2) BYETA=BYETA/2**(IB-1) 
        BPHI=PARU(2) 
        IF(IM.NE.1) BPHI=BPHI/2**(IB-1) 
        IF(IM.LE.2) BNAVE=FAC*NMUFM/FLOAT(2**(IB-1)) 
        IF(IM.EQ.3) BNAVE=FAC*NMUFM/FLOAT(4**(IB-1)) 
        DO 490 IP=1,4 
        FMOMA(IP)=FAC*FM1FM(IM,IB,IP) 
        FMOMS(IP)=SQRT(MAX(0.,FAC*(FAC*FM2FM(IM,IB,IP)-FMOMA(IP)**2))) 
  490   CONTINUE 
        WRITE(MSTU(11),5600) BYETA,BPHI,BNAVE,(FMOMA(IP),FMOMS(IP), 
     &  IP=1,4) 
  500   CONTINUE 
  510   CONTINUE 
 
C...Copy statistics on factorial moments into /LUJETS/. 
      ELSEIF(MTABU.EQ.33) THEN 
        FAC=1./MAX(1,NEVFM) 
        DO 540 IM=1,3 
        DO 530 IB=1,10 
        I=10*(IM-1)+IB 
        K(I,1)=32 
        K(I,2)=99 
        K(I,3)=1 
        IF(IM.NE.2) K(I,3)=2**(IB-1) 
        K(I,4)=1 
        IF(IM.NE.1) K(I,4)=2**(IB-1) 
        K(I,5)=0 
        P(I,1)=2.*PARU(57)/K(I,3) 
        V(I,1)=PARU(2)/K(I,4) 
        DO 520 IP=1,4 
        P(I,IP+1)=FAC*FM1FM(IM,IB,IP) 
        V(I,IP+1)=SQRT(MAX(0.,FAC*(FAC*FM2FM(IM,IB,IP)-P(I,IP+1)**2))) 
  520   CONTINUE 
  530   CONTINUE 
  540   CONTINUE 
        N=30 
        DO 550 J=1,5 
        K(N+1,J)=0 
        P(N+1,J)=0. 
        V(N+1,J)=0. 
  550   CONTINUE 
        K(N+1,1)=32 
        K(N+1,2)=99 
        K(N+1,5)=NEVFM 
        MSTU(3)=1 
 
C...Reset statistics on Energy-Energy Correlation. 
      ELSEIF(MTABU.EQ.40) THEN 
        NEVEE=0 
        DO 560 J=1,25 
        FE1EC(J)=0. 
        FE2EC(J)=0. 
        FE1EC(51-J)=0. 
        FE2EC(51-J)=0. 
        FE1EA(J)=0. 
        FE2EA(J)=0. 
  560   CONTINUE 
 
C...Find particles to include, with proper assumed mass. 
      ELSEIF(MTABU.EQ.41) THEN 
        NEVEE=NEVEE+1 
        NLOW=N+MSTU(3) 
        NUPP=NLOW 
        ECM=0. 
        DO 570 I=1,N 
        IF(K(I,1).LE.0.OR.K(I,1).GT.10) GOTO 570 
        IF(MSTU(41).GE.2) THEN 
          KC=LUCOMP(K(I,2)) 
          IF(KC.EQ.0.OR.KC.EQ.12.OR.KC.EQ.14.OR.KC.EQ.16.OR. 
     &    KC.EQ.18) GOTO 570 
          IF(MSTU(41).GE.3.AND.KCHG(KC,2).EQ.0.AND.LUCHGE(K(I,2)).EQ.0) 
     &    GOTO 570 
        ENDIF 
        PMR=0. 
        IF(MSTU(42).EQ.1.AND.K(I,2).NE.22) PMR=ULMASS(211) 
        IF(MSTU(42).GE.2) PMR=P(I,5) 
        IF(NUPP.GT.MSTU(4)-5-MSTU(32)) THEN 
          CALL LUERRM(11,'(LUTABU:) no more memory left in LUJETS') 
          RETURN 
        ENDIF 
        NUPP=NUPP+1 
        P(NUPP,1)=P(I,1) 
        P(NUPP,2)=P(I,2) 
        P(NUPP,3)=P(I,3) 
        P(NUPP,4)=SQRT(PMR**2+P(I,1)**2+P(I,2)**2+P(I,3)**2) 
        P(NUPP,5)=MAX(1E-10,SQRT(P(I,1)**2+P(I,2)**2+P(I,3)**2)) 
        ECM=ECM+P(NUPP,4) 
  570   CONTINUE 
        IF(NUPP.EQ.NLOW) RETURN 
 
C...Analyze Energy-Energy Correlation in event. 
        FAC=(2./ECM**2)*50./PARU(1) 
        DO 580 J=1,50 
        FEVEE(J)=0. 
  580   CONTINUE 
        DO 600 I1=NLOW+2,NUPP 
        DO 590 I2=NLOW+1,I1-1 
        CTHE=(P(I1,1)*P(I2,1)+P(I1,2)*P(I2,2)+P(I1,3)*P(I2,3))/ 
     &  (P(I1,5)*P(I2,5)) 
        THE=ACOS(MAX(-1.,MIN(1.,CTHE))) 
        ITHE=MAX(1,MIN(50,1+INT(50.*THE/PARU(1)))) 
        FEVEE(ITHE)=FEVEE(ITHE)+FAC*P(I1,4)*P(I2,4) 
  590   CONTINUE 
  600   CONTINUE 
        DO 610 J=1,25 
        FE1EC(J)=FE1EC(J)+FEVEE(J) 
        FE2EC(J)=FE2EC(J)+FEVEE(J)**2 
        FE1EC(51-J)=FE1EC(51-J)+FEVEE(51-J) 
        FE2EC(51-J)=FE2EC(51-J)+FEVEE(51-J)**2 
        FE1EA(J)=FE1EA(J)+(FEVEE(51-J)-FEVEE(J)) 
        FE2EA(J)=FE2EA(J)+(FEVEE(51-J)-FEVEE(J))**2 
  610   CONTINUE 
        MSTU(62)=NUPP-NLOW 
 
C...Write statistics on Energy-Energy Correlation. 
      ELSEIF(MTABU.EQ.42) THEN 
        FAC=1./MAX(1,NEVEE) 
        WRITE(MSTU(11),5700) NEVEE 
        DO 620 J=1,25 
        FEEC1=FAC*FE1EC(J) 
        FEES1=SQRT(MAX(0.,FAC*(FAC*FE2EC(J)-FEEC1**2))) 
        FEEC2=FAC*FE1EC(51-J) 
        FEES2=SQRT(MAX(0.,FAC*(FAC*FE2EC(51-J)-FEEC2**2))) 
        FEECA=FAC*FE1EA(J) 
        FEESA=SQRT(MAX(0.,FAC*(FAC*FE2EA(J)-FEECA**2))) 
        WRITE(MSTU(11),5800) 3.6*(J-1),3.6*J,FEEC1,FEES1,FEEC2,FEES2, 
     &  FEECA,FEESA 
  620   CONTINUE 
 
C...Copy statistics on Energy-Energy Correlation into /LUJETS/. 
      ELSEIF(MTABU.EQ.43) THEN 
        FAC=1./MAX(1,NEVEE) 
        DO 630 I=1,25 
        K(I,1)=32 
        K(I,2)=99 
        K(I,3)=0 
        K(I,4)=0 
        K(I,5)=0 
        P(I,1)=FAC*FE1EC(I) 
        V(I,1)=SQRT(MAX(0.,FAC*(FAC*FE2EC(I)-P(I,1)**2))) 
        P(I,2)=FAC*FE1EC(51-I) 
        V(I,2)=SQRT(MAX(0.,FAC*(FAC*FE2EC(51-I)-P(I,2)**2))) 
        P(I,3)=FAC*FE1EA(I) 
        V(I,3)=SQRT(MAX(0.,FAC*(FAC*FE2EA(I)-P(I,3)**2))) 
        P(I,4)=PARU(1)*(I-1)/50. 
        P(I,5)=PARU(1)*I/50. 
        V(I,4)=3.6*(I-1) 
        V(I,5)=3.6*I 
  630   CONTINUE 
        N=25 
        DO 640 J=1,5 
        K(N+1,J)=0 
        P(N+1,J)=0. 
        V(N+1,J)=0. 
  640   CONTINUE 
        K(N+1,1)=32 
        K(N+1,2)=99 
        K(N+1,5)=NEVEE 
        MSTU(3)=1 
 
C...Reset statistics on decay channels. 
      ELSEIF(MTABU.EQ.50) THEN 
        NEVDC=0 
        NKFDC=0 
        NREDC=0 
 
C...Identify and order flavour content of final state. 
      ELSEIF(MTABU.EQ.51) THEN 
        NEVDC=NEVDC+1 
        NDS=0 
        DO 670 I=1,N 
        IF(K(I,1).LE.0.OR.K(I,1).GE.6) GOTO 670 
        NDS=NDS+1 
        IF(NDS.GT.8) THEN 
          NREDC=NREDC+1 
          RETURN 
        ENDIF 
        KFM=2*IABS(K(I,2)) 
        IF(K(I,2).LT.0) KFM=KFM-1 
        DO 650 IDS=NDS-1,1,-1 
        IIN=IDS+1 
        IF(KFM.LT.KFDM(IDS)) GOTO 660 
        KFDM(IDS+1)=KFDM(IDS) 
  650   CONTINUE 
        IIN=1 
  660   KFDM(IIN)=KFM 
  670   CONTINUE 
 
C...Find whether old or new final state. 
        DO 690 IDC=1,NKFDC 
        IF(NDS.LT.KFDC(IDC,0)) THEN 
          IKFDC=IDC 
          GOTO 700 
        ELSEIF(NDS.EQ.KFDC(IDC,0)) THEN 
          DO 680 I=1,NDS 
          IF(KFDM(I).LT.KFDC(IDC,I)) THEN 
            IKFDC=IDC 
            GOTO 700 
          ELSEIF(KFDM(I).GT.KFDC(IDC,I)) THEN 
            GOTO 690 
          ENDIF 
  680     CONTINUE 
          IKFDC=-IDC 
          GOTO 700 
        ENDIF 
  690   CONTINUE 
        IKFDC=NKFDC+1 
  700   IF(IKFDC.LT.0) THEN 
          IKFDC=-IKFDC 
        ELSEIF(NKFDC.GE.200) THEN 
          NREDC=NREDC+1 
          RETURN 
        ELSE 
          DO 720 IDC=NKFDC,IKFDC,-1 
          NPDC(IDC+1)=NPDC(IDC) 
          DO 710 I=0,8 
          KFDC(IDC+1,I)=KFDC(IDC,I) 
  710     CONTINUE 
  720     CONTINUE 
          NKFDC=NKFDC+1 
          KFDC(IKFDC,0)=NDS 
          DO 730 I=1,NDS 
          KFDC(IKFDC,I)=KFDM(I) 
  730     CONTINUE 
          NPDC(IKFDC)=0 
        ENDIF 
        NPDC(IKFDC)=NPDC(IKFDC)+1 
 
C...Write statistics on decay channels. 
      ELSEIF(MTABU.EQ.52) THEN 
        FAC=1./MAX(1,NEVDC) 
        WRITE(MSTU(11),5900) NEVDC 
        DO 750 IDC=1,NKFDC 
        DO 740 I=1,KFDC(IDC,0) 
        KFM=KFDC(IDC,I) 
        KF=(KFM+1)/2 
        IF(2*KF.NE.KFM) KF=-KF 
        CALL LUNAME(KF,CHAU) 
        CHDC(I)=CHAU(1:12) 
        IF(CHAU(13:13).NE.' ') CHDC(I)(12:12)='?' 
  740   CONTINUE 
        WRITE(MSTU(11),6000) FAC*NPDC(IDC),(CHDC(I),I=1,KFDC(IDC,0)) 
  750   CONTINUE 
        IF(NREDC.NE.0) WRITE(MSTU(11),6100) FAC*NREDC 
 
C...Copy statistics on decay channels into /LUJETS/. 
      ELSEIF(MTABU.EQ.53) THEN 
        FAC=1./MAX(1,NEVDC) 
        DO 780 IDC=1,NKFDC 
        K(IDC,1)=32 
        K(IDC,2)=99 
        K(IDC,3)=0 
        K(IDC,4)=0 
        K(IDC,5)=KFDC(IDC,0) 
        DO 760 J=1,5 
        P(IDC,J)=0. 
        V(IDC,J)=0. 
  760   CONTINUE 
        DO 770 I=1,KFDC(IDC,0) 
        KFM=KFDC(IDC,I) 
        KF=(KFM+1)/2 
        IF(2*KF.NE.KFM) KF=-KF 
        IF(I.LE.5) P(IDC,I)=KF 
        IF(I.GE.6) V(IDC,I-5)=KF 
  770   CONTINUE 
        V(IDC,5)=FAC*NPDC(IDC) 
  780   CONTINUE 
        N=NKFDC 
        DO 790 J=1,5 
        K(N+1,J)=0 
        P(N+1,J)=0. 
        V(N+1,J)=0. 
  790   CONTINUE 
        K(N+1,1)=32 
        K(N+1,2)=99 
        K(N+1,5)=NEVDC 
        V(N+1,5)=FAC*NREDC 
        MSTU(3)=1 
      ENDIF 
 
C...Format statements for output on unit MSTU(11) (default 6). 
 5000 FORMAT(///20X,'Event statistics - initial state'/ 
     &20X,'based on an analysis of ',I6,' events'// 
     &3X,'Main flavours after',8X,'Fraction',4X,'Subfractions ', 
     &'according to fragmenting system multiplicity'/ 
     &4X,'hard interaction',24X,'1',7X,'2',7X,'3',7X,'4',7X,'5', 
     &6X,'6-7',5X,'8-10',3X,'11-15',3X,'16-25',4X,'>25'/) 
 5100 FORMAT(3X,A12,1X,A12,F10.5,1X,10F8.4) 
 5200 FORMAT(///20X,'Event statistics - final state'/ 
     &20X,'based on an analysis of ',I7,' events'// 
     &5X,'Mean primary multiplicity =',F10.4/ 
     &5X,'Mean final   multiplicity =',F10.4/ 
     &5X,'Mean charged multiplicity =',F10.4// 
     &5X,'Number of particles produced per event (directly and via ', 
     &'decays/branchings)'/ 
     &5X,'KF    Particle/jet  MDCY',10X,'Particles',13X,'Antiparticles', 
     &8X,'Total'/35X,'prim        seco        prim        seco'/) 
 5300 FORMAT(1X,I6,4X,A16,I2,5(1X,F11.6)) 
 5400 FORMAT(///20X,'Factorial moments analysis of multiplicity'/ 
     &20X,'based on an analysis of ',I6,' events'// 
     &3X,'delta-',A3,' delta-phi     <n>/bin',10X,'<F2>',18X,'<F3>', 
     &18X,'<F4>',18X,'<F5>'/35X,4('     value     error  ')) 
 5500 FORMAT(10X) 
 5600 FORMAT(2X,2F10.4,F12.4,4(F12.4,F10.4)) 
 5700 FORMAT(///20X,'Energy-Energy Correlation and Asymmetry'/ 
     &20X,'based on an analysis of ',I6,' events'// 
     &2X,'theta range',8X,'EEC(theta)',8X,'EEC(180-theta)',7X, 
     &'EECA(theta)'/2X,'in degrees ',3('      value    error')/) 
 5800 FORMAT(2X,F4.1,' - ',F4.1,3(F11.4,F9.4)) 
 5900 FORMAT(///20X,'Decay channel analysis - final state'/ 
     &20X,'based on an analysis of ',I6,' events'// 
     &2X,'Probability',10X,'Complete final state'/) 
 6000 FORMAT(2X,F9.5,5X,8(A12,1X)) 
 6100 FORMAT(2X,F9.5,5X,'into other channels (more than 8 particles ', 
     &'or table overflow)') 
 
      RETURN 
      END 
