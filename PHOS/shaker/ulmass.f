*CMZ :          17/07/98  15.44.33  by  Federico Carminati
*-- Author :
C*********************************************************************

      FUNCTION ULMASS(KF)

C...Purpose: to give the mass of a particle/parton.
*KEEP,LUDAT1.
      COMMON /LUDAT1/ MSTU(200),PARU(200),MSTJ(200),PARJ(200)
      SAVE /LUDAT1/
*KEEP,LUDAT2.
      COMMON /LUDAT2/ KCHG(500,3),PMAS(500,4),PARF(2000),VCKM(4,4)
      SAVE /LUDAT2/
*KEND.

C...Reset variables. Compressed code.
      ULMASS=0.
      KFA=IABS(KF)
      KC=LUCOMP(KF)
      IF(KC.EQ.0) RETURN
      PARF(106)=PMAS(6,1)
      PARF(107)=PMAS(7,1)
      PARF(108)=PMAS(8,1)

C...Guarantee use of constituent masses for internal checks.
      IF((MSTJ(93).EQ.1.OR.MSTJ(93).EQ.2).AND.KFA.LE.10) THEN
        ULMASS=PARF(100+KFA)
        IF(MSTJ(93).EQ.2) ULMASS=MAX(0.,ULMASS-PARF(121))

C...Masses that can be read directly off table.
      ELSEIF(KFA.LE.100.OR.KC.LE.80.OR.KC.GT.100) THEN
        ULMASS=PMAS(KC,1)

C...Find constituent partons and their masses.
      ELSE
        KFLA=MOD(KFA/1000,10)
        KFLB=MOD(KFA/100,10)
        KFLC=MOD(KFA/10,10)
        KFLS=MOD(KFA,10)
        KFLR=MOD(KFA/10000,10)
        PMA=PARF(100+KFLA)
        PMB=PARF(100+KFLB)
        PMC=PARF(100+KFLC)

C...Construct masses for various meson, diquark and baryon cases.
        IF(KFLA.EQ.0.AND.KFLR.EQ.0.AND.KFLS.LE.3) THEN
          IF(KFLS.EQ.1) PMSPL=-3./(PMB*PMC)
          IF(KFLS.GE.3) PMSPL=1./(PMB*PMC)
          ULMASS=PARF(111)+PMB+PMC+PARF(113)*PARF(101)**2*PMSPL
        ELSEIF(KFLA.EQ.0) THEN
          KMUL=2
          IF(KFLS.EQ.1) KMUL=3
          IF(KFLR.EQ.2) KMUL=4
          IF(KFLS.EQ.5) KMUL=5
          ULMASS=PARF(113+KMUL)+PMB+PMC
        ELSEIF(KFLC.EQ.0) THEN
          IF(KFLS.EQ.1) PMSPL=-3./(PMA*PMB)
          IF(KFLS.EQ.3) PMSPL=1./(PMA*PMB)
          ULMASS=2.*PARF(112)/3.+PMA+PMB+PARF(114)*PARF(101)**2*PMSPL
          IF(MSTJ(93).EQ.1) ULMASS=PMA+PMB
          IF(MSTJ(93).EQ.2) ULMASS=MAX(0.,ULMASS-PARF(122)-
     &    2.*PARF(112)/3.)
        ELSE
          IF(KFLS.EQ.2.AND.KFLA.EQ.KFLB) THEN
            PMSPL=1./(PMA*PMB)-2./(PMA*PMC)-2./(PMB*PMC)
          ELSEIF(KFLS.EQ.2.AND.KFLB.GE.KFLC) THEN
            PMSPL=-2./(PMA*PMB)-2./(PMA*PMC)+1./(PMB*PMC)
          ELSEIF(KFLS.EQ.2) THEN
            PMSPL=-3./(PMB*PMC)
          ELSE
            PMSPL=1./(PMA*PMB)+1./(PMA*PMC)+1./(PMB*PMC)
          ENDIF
          ULMASS=PARF(112)+PMA+PMB+PMC+PARF(114)*PARF(101)**2*PMSPL
        ENDIF
      ENDIF

C...Optional mass broadening according to truncated Breit-Wigner
C...(either in m or in m^2).
      IF(MSTJ(24).GE.1.AND.PMAS(KC,2).GT.1E-4) THEN
        IF(MSTJ(24).EQ.1.OR.(MSTJ(24).EQ.2.AND.KFA.GT.100)) THEN
          ULMASS=ULMASS+0.5*PMAS(KC,2)*TAN((2.*RLU(0)-1.)*
     &    ATAN(2.*PMAS(KC,3)/PMAS(KC,2)))
        ELSE
          PM0=ULMASS
          PMLOW=ATAN((MAX(0.,PM0-PMAS(KC,3))**2-PM0**2)/
     &    (PM0*PMAS(KC,2)))
          PMUPP=ATAN(((PM0+PMAS(KC,3))**2-PM0**2)/(PM0*PMAS(KC,2)))
          ULMASS=SQRT(MAX(0.,PM0**2+PM0*PMAS(KC,2)*TAN(PMLOW+
     &    (PMUPP-PMLOW)*RLU(0))))
        ENDIF
      ENDIF
      MSTJ(93)=0

      RETURN
      END
