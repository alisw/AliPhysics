*CMZ :          17/07/98  15.44.32  by  Federico Carminati
*-- Author :
C*********************************************************************

      SUBROUTINE LUKFDI(KFL1,KFL2,KFL3,KF)

C...Purpose: to generate a new flavour pair and combine off a hadron.
*KEEP,LUDAT1.
      COMMON /LUDAT1/ MSTU(200),PARU(200),MSTJ(200),PARJ(200)
      SAVE /LUDAT1/
*KEEP,LUDAT2.
      COMMON /LUDAT2/ KCHG(500,3),PMAS(500,4),PARF(2000),VCKM(4,4)
      SAVE /LUDAT2/
*KEND.

C...Default flavour values. Input consistency checks.
      KF1A=IABS(KFL1)
      KF2A=IABS(KFL2)
      KFL3=0
      KF=0
      IF(KF1A.EQ.0) RETURN
      IF(KF2A.NE.0) THEN
        IF(KF1A.LE.10.AND.KF2A.LE.10.AND.KFL1*KFL2.GT.0) RETURN
        IF(KF1A.GT.10.AND.KF2A.GT.10) RETURN
        IF((KF1A.GT.10.OR.KF2A.GT.10).AND.KFL1*KFL2.LT.0) RETURN
      ENDIF

C...Check if tabulated flavour probabilities are to be used.
      IF(MSTJ(15).EQ.1) THEN
        KTAB1=-1
        IF(KF1A.GE.1.AND.KF1A.LE.6) KTAB1=KF1A
        KFL1A=MOD(KF1A/1000,10)
        KFL1B=MOD(KF1A/100,10)
        KFL1S=MOD(KF1A,10)
        IF(KFL1A.GE.1.AND.KFL1A.LE.4.AND.KFL1B.GE.1.AND.KFL1B.LE.4)
     &  KTAB1=6+KFL1A*(KFL1A-2)+2*KFL1B+(KFL1S-1)/2
        IF(KFL1A.GE.1.AND.KFL1A.LE.4.AND.KFL1A.EQ.KFL1B) KTAB1=KTAB1-1
        IF(KF1A.GE.1.AND.KF1A.LE.6) KFL1A=KF1A
        KTAB2=0
        IF(KF2A.NE.0) THEN
          KTAB2=-1
          IF(KF2A.GE.1.AND.KF2A.LE.6) KTAB2=KF2A
          KFL2A=MOD(KF2A/1000,10)
          KFL2B=MOD(KF2A/100,10)
          KFL2S=MOD(KF2A,10)
          IF(KFL2A.GE.1.AND.KFL2A.LE.4.AND.KFL2B.GE.1.AND.KFL2B.LE.4)
     &    KTAB2=6+KFL2A*(KFL2A-2)+2*KFL2B+(KFL2S-1)/2
          IF(KFL2A.GE.1.AND.KFL2A.LE.4.AND.KFL2A.EQ.KFL2B) KTAB2=KTAB2-1
        ENDIF
        IF(KTAB1.GE.0.AND.KTAB2.GE.0) GOTO 140
      ENDIF

C...Parameters and breaking diquark parameter combinations.
  100 PAR2=PARJ(2)
      PAR3=PARJ(3)
      PAR4=3.*PARJ(4)
      IF(MSTJ(12).GE.2) THEN
        PAR3M=SQRT(PARJ(3))
        PAR4M=1./(3.*SQRT(PARJ(4)))
        PARDM=PARJ(7)/(PARJ(7)+PAR3M*PARJ(6))
        PARS0=PARJ(5)*(2.+(1.+PAR2*PAR3M*PARJ(7))*(1.+PAR4M))
        PARS1=PARJ(7)*PARS0/(2.*PAR3M)+PARJ(5)*(PARJ(6)*(1.+PAR4M)+
     &  PAR2*PAR3M*PARJ(6)*PARJ(7))
        PARS2=PARJ(5)*2.*PARJ(6)*PARJ(7)*(PAR2*PARJ(7)+(1.+PAR4M)/PAR3M)
        PARSM=MAX(PARS0,PARS1,PARS2)
        PAR4=PAR4*(1.+PARSM)/(1.+PARSM/(3.*PAR4M))
      ENDIF

C...Choice of whether to generate meson or baryon.
      MBARY=0
      KFDA=0
      IF(KF1A.LE.10) THEN
        IF(KF2A.EQ.0.AND.MSTJ(12).GE.1.AND.(1.+PARJ(1))*RLU(0).GT.1.)
     &  MBARY=1
        IF(KF2A.GT.10) MBARY=2
        IF(KF2A.GT.10.AND.KF2A.LE.10000) KFDA=KF2A
      ELSE
        MBARY=2
        IF(KF1A.LE.10000) KFDA=KF1A
      ENDIF

C...Possibility of process diquark -> meson + new diquark.
      IF(KFDA.NE.0.AND.MSTJ(12).GE.2) THEN
        KFLDA=MOD(KFDA/1000,10)
        KFLDB=MOD(KFDA/100,10)
        KFLDS=MOD(KFDA,10)
        WTDQ=PARS0
        IF(MAX(KFLDA,KFLDB).EQ.3) WTDQ=PARS1
        IF(MIN(KFLDA,KFLDB).EQ.3) WTDQ=PARS2
        IF(KFLDS.EQ.1) WTDQ=WTDQ/(3.*PAR4M)
        IF((1.+WTDQ)*RLU(0).GT.1.) MBARY=-1
        IF(MBARY.EQ.-1.AND.KF2A.NE.0) RETURN
      ENDIF

C...Flavour for meson, possibly with new flavour.
      IF(MBARY.LE.0) THEN
        KFS=ISIGN(1,KFL1)
        IF(MBARY.EQ.0) THEN
          IF(KF2A.EQ.0) KFL3=ISIGN(1+INT((2.+PAR2)*RLU(0)),-KFL1)
          KFLA=MAX(KF1A,KF2A+IABS(KFL3))
          KFLB=MIN(KF1A,KF2A+IABS(KFL3))
          IF(KFLA.NE.KF1A) KFS=-KFS

C...Splitting of diquark into meson plus new diquark.
        ELSE
          KFL1A=MOD(KF1A/1000,10)
          KFL1B=MOD(KF1A/100,10)
  110     KFL1D=KFL1A+INT(RLU(0)+0.5)*(KFL1B-KFL1A)
          KFL1E=KFL1A+KFL1B-KFL1D
          IF((KFL1D.EQ.3.AND.RLU(0).GT.PARDM).OR.(KFL1E.EQ.3.AND.
     &    RLU(0).LT.PARDM)) THEN
            KFL1D=KFL1A+KFL1B-KFL1D
            KFL1E=KFL1A+KFL1B-KFL1E
          ENDIF
          KFL3A=1+INT((2.+PAR2*PAR3M*PARJ(7))*RLU(0))
          IF((KFL1E.NE.KFL3A.AND.RLU(0).GT.(1.+PAR4M)/MAX(2.,1.+PAR4M)).
     &    OR.(KFL1E.EQ.KFL3A.AND.RLU(0).GT.2./MAX(2.,1.+PAR4M)))
     &    GOTO 110
          KFLDS=3
          IF(KFL1E.NE.KFL3A) KFLDS=2*INT(RLU(0)+1./(1.+PAR4M))+1
          KFL3=ISIGN(10000+1000*MAX(KFL1E,KFL3A)+100*MIN(KFL1E,KFL3A)+
     &    KFLDS,-KFL1)
          KFLA=MAX(KFL1D,KFL3A)
          KFLB=MIN(KFL1D,KFL3A)
          IF(KFLA.NE.KFL1D) KFS=-KFS
        ENDIF

C...Form meson, with spin and flavour mixing for diagonal states.
        IF(KFLA.LE.2) KMUL=INT(PARJ(11)+RLU(0))
        IF(KFLA.EQ.3) KMUL=INT(PARJ(12)+RLU(0))
        IF(KFLA.GE.4) KMUL=INT(PARJ(13)+RLU(0))
        IF(KMUL.EQ.0.AND.PARJ(14).GT.0.) THEN
          IF(RLU(0).LT.PARJ(14)) KMUL=2
        ELSEIF(KMUL.EQ.1.AND.PARJ(15)+PARJ(16)+PARJ(17).GT.0.) THEN
          RMUL=RLU(0)
          IF(RMUL.LT.PARJ(15)) KMUL=3
          IF(KMUL.EQ.1.AND.RMUL.LT.PARJ(15)+PARJ(16)) KMUL=4
          IF(KMUL.EQ.1.AND.RMUL.LT.PARJ(15)+PARJ(16)+PARJ(17)) KMUL=5
        ENDIF
        KFLS=3
        IF(KMUL.EQ.0.OR.KMUL.EQ.3) KFLS=1
        IF(KMUL.EQ.5) KFLS=5
        IF(KFLA.NE.KFLB) THEN
          KF=(100*KFLA+10*KFLB+KFLS)*KFS*(-1)**KFLA
        ELSE
          RMIX=RLU(0)
          IMIX=2*KFLA+10*KMUL
          IF(KFLA.LE.3) KF=110*(1+INT(RMIX+PARF(IMIX-1))+
     &    INT(RMIX+PARF(IMIX)))+KFLS
          IF(KFLA.GE.4) KF=110*KFLA+KFLS
        ENDIF
        IF(KMUL.EQ.2.OR.KMUL.EQ.3) KF=KF+ISIGN(10000,KF)
        IF(KMUL.EQ.4) KF=KF+ISIGN(20000,KF)

C...Generate diquark flavour.
      ELSE
  120   IF(KF1A.LE.10.AND.KF2A.EQ.0) THEN
          KFLA=KF1A
  130     KFLB=1+INT((2.+PAR2*PAR3)*RLU(0))
          KFLC=1+INT((2.+PAR2*PAR3)*RLU(0))
          KFLDS=1
          IF(KFLB.GE.KFLC) KFLDS=3
          IF(KFLDS.EQ.1.AND.PAR4*RLU(0).GT.1.) GOTO 130
          IF(KFLDS.EQ.3.AND.PAR4.LT.RLU(0)) GOTO 130
          KFL3=ISIGN(1000*MAX(KFLB,KFLC)+100*MIN(KFLB,KFLC)+KFLDS,KFL1)

C...Take diquark flavour from input.
        ELSEIF(KF1A.LE.10) THEN
          KFLA=KF1A
          KFLB=MOD(KF2A/1000,10)
          KFLC=MOD(KF2A/100,10)
          KFLDS=MOD(KF2A,10)

C...Generate (or take from input) quark to go with diquark.
        ELSE
          IF(KF2A.EQ.0) KFL3=ISIGN(1+INT((2.+PAR2)*RLU(0)),KFL1)
          KFLA=KF2A+IABS(KFL3)
          KFLB=MOD(KF1A/1000,10)
          KFLC=MOD(KF1A/100,10)
          KFLDS=MOD(KF1A,10)
        ENDIF

C...SU(6) factors for formation of baryon. Try again if fails.
        KBARY=KFLDS
        IF(KFLDS.EQ.3.AND.KFLB.NE.KFLC) KBARY=5
        IF(KFLA.NE.KFLB.AND.KFLA.NE.KFLC) KBARY=KBARY+1
        WT=PARF(60+KBARY)+PARJ(18)*PARF(70+KBARY)
        IF(MBARY.EQ.1.AND.MSTJ(12).GE.2) THEN
          WTDQ=PARS0
          IF(MAX(KFLB,KFLC).EQ.3) WTDQ=PARS1
          IF(MIN(KFLB,KFLC).EQ.3) WTDQ=PARS2
          IF(KFLDS.EQ.1) WTDQ=WTDQ/(3.*PAR4M)
          IF(KFLDS.EQ.1) WT=WT*(1.+WTDQ)/(1.+PARSM/(3.*PAR4M))
          IF(KFLDS.EQ.3) WT=WT*(1.+WTDQ)/(1.+PARSM)
        ENDIF
        IF(KF2A.EQ.0.AND.WT.LT.RLU(0)) GOTO 120

C...Form baryon. Distinguish Lambda- and Sigmalike baryons.
        KFLD=MAX(KFLA,KFLB,KFLC)
        KFLF=MIN(KFLA,KFLB,KFLC)
        KFLE=KFLA+KFLB+KFLC-KFLD-KFLF
        KFLS=2
        IF((PARF(60+KBARY)+PARJ(18)*PARF(70+KBARY))*RLU(0).GT.
     &  PARF(60+KBARY)) KFLS=4
        KFLL=0
        IF(KFLS.EQ.2.AND.KFLD.GT.KFLE.AND.KFLE.GT.KFLF) THEN
          IF(KFLDS.EQ.1.AND.KFLA.EQ.KFLD) KFLL=1
          IF(KFLDS.EQ.1.AND.KFLA.NE.KFLD) KFLL=INT(0.25+RLU(0))
          IF(KFLDS.EQ.3.AND.KFLA.NE.KFLD) KFLL=INT(0.75+RLU(0))
        ENDIF
        IF(KFLL.EQ.0) KF=ISIGN(1000*KFLD+100*KFLE+10*KFLF+KFLS,KFL1)
        IF(KFLL.EQ.1) KF=ISIGN(1000*KFLD+100*KFLF+10*KFLE+KFLS,KFL1)
      ENDIF
      RETURN

C...Use tabulated probabilities to select new flavour and hadron.
  140 IF(KTAB2.EQ.0.AND.MSTJ(12).LE.0) THEN
        KT3L=1
        KT3U=6
      ELSEIF(KTAB2.EQ.0.AND.KTAB1.GE.7.AND.MSTJ(12).LE.1) THEN
        KT3L=1
        KT3U=6
      ELSEIF(KTAB2.EQ.0) THEN
        KT3L=1
        KT3U=22
      ELSE
        KT3L=KTAB2
        KT3U=KTAB2
      ENDIF
      RFL=0.
      DO 150 KTS=0,2
      DO 150 KT3=KT3L,KT3U
      RFL=RFL+PARF(120+80*KTAB1+25*KTS+KT3)
  150 CONTINUE
      RFL=RLU(0)*RFL
      DO 160 KTS=0,2
      KTABS=KTS
      DO 160 KT3=KT3L,KT3U
      KTAB3=KT3
      RFL=RFL-PARF(120+80*KTAB1+25*KTS+KT3)
  160 IF(RFL.LE.0.) GOTO 170
  170 CONTINUE

C...Reconstruct flavour of produced quark/diquark.
      IF(KTAB3.LE.6) THEN
        KFL3A=KTAB3
        KFL3B=0
        KFL3=ISIGN(KFL3A,KFL1*(2*KTAB1-13))
      ELSE
        KFL3A=1
        IF(KTAB3.GE.8) KFL3A=2
        IF(KTAB3.GE.11) KFL3A=3
        IF(KTAB3.GE.16) KFL3A=4
        KFL3B=(KTAB3-6-KFL3A*(KFL3A-2))/2
        KFL3=1000*KFL3A+100*KFL3B+1
        IF(KFL3A.EQ.KFL3B.OR.KTAB3.NE.6+KFL3A*(KFL3A-2)+2*KFL3B) KFL3=
     &  KFL3+2
        KFL3=ISIGN(KFL3,KFL1*(13-2*KTAB1))
      ENDIF

C...Reconstruct meson code.
      IF(KFL3A.EQ.KFL1A.AND.KFL3B.EQ.KFL1B.AND.(KFL3A.LE.3.OR.
     &KFL3B.NE.0)) THEN
        RFL=RLU(0)*(PARF(143+80*KTAB1+25*KTABS)+PARF(144+80*KTAB1+
     &  25*KTABS)+PARF(145+80*KTAB1+25*KTABS))
        KF=110+2*KTABS+1
        IF(RFL.GT.PARF(143+80*KTAB1+25*KTABS)) KF=220+2*KTABS+1
        IF(RFL.GT.PARF(143+80*KTAB1+25*KTABS)+PARF(144+80*KTAB1+
     &  25*KTABS)) KF=330+2*KTABS+1
      ELSEIF(KTAB1.LE.6.AND.KTAB3.LE.6) THEN
        KFLA=MAX(KTAB1,KTAB3)
        KFLB=MIN(KTAB1,KTAB3)
        KFS=ISIGN(1,KFL1)
        IF(KFLA.NE.KF1A) KFS=-KFS
        KF=(100*KFLA+10*KFLB+2*KTABS+1)*KFS*(-1)**KFLA
      ELSEIF(KTAB1.GE.7.AND.KTAB3.GE.7) THEN
        KFS=ISIGN(1,KFL1)
        IF(KFL1A.EQ.KFL3A) THEN
          KFLA=MAX(KFL1B,KFL3B)
          KFLB=MIN(KFL1B,KFL3B)
          IF(KFLA.NE.KFL1B) KFS=-KFS
        ELSEIF(KFL1A.EQ.KFL3B) THEN
          KFLA=KFL3A
          KFLB=KFL1B
          KFS=-KFS
        ELSEIF(KFL1B.EQ.KFL3A) THEN
          KFLA=KFL1A
          KFLB=KFL3B
        ELSEIF(KFL1B.EQ.KFL3B) THEN
          KFLA=MAX(KFL1A,KFL3A)
          KFLB=MIN(KFL1A,KFL3A)
          IF(KFLA.NE.KFL1A) KFS=-KFS
        ELSE
          CALL LUERRM(2,'(LUKFDI:) no matching flavours for qq -> qq')
          GOTO 100
        ENDIF
        KF=(100*KFLA+10*KFLB+2*KTABS+1)*KFS*(-1)**KFLA

C...Reconstruct baryon code.
      ELSE
        IF(KTAB1.GE.7) THEN
          KFLA=KFL3A
          KFLB=KFL1A
          KFLC=KFL1B
        ELSE
          KFLA=KFL1A
          KFLB=KFL3A
          KFLC=KFL3B
        ENDIF
        KFLD=MAX(KFLA,KFLB,KFLC)
        KFLF=MIN(KFLA,KFLB,KFLC)
        KFLE=KFLA+KFLB+KFLC-KFLD-KFLF
        IF(KTABS.EQ.0) KF=ISIGN(1000*KFLD+100*KFLF+10*KFLE+2,KFL1)
        IF(KTABS.GE.1) KF=ISIGN(1000*KFLD+100*KFLE+10*KFLF+2*KTABS,KFL1)
      ENDIF

C...Check that constructed flavour code is an allowed one.
      IF(KFL2.NE.0) KFL3=0
      KC=LUCOMP(KF)
      IF(KC.EQ.0) THEN
        CALL LUERRM(2,'(LUKFDI:) user-defined flavour probabilities '//
     &  'failed')
        GOTO 100
      ENDIF

      RETURN
      END
