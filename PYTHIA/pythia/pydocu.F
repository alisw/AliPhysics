 
C*********************************************************************
 
      SUBROUTINE PYDOCU
 
C...Handles the decumentation of the process in MSTI and PARI,
C...and also computes cross-sections based on accumulated statistics.
      COMMON/LUJETS/N,K(4000,5),P(4000,5),V(4000,5)
      COMMON/LUDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
      COMMON/PYPARS/MSTP(200),PARP(200),MSTI(200),PARI(200)
      COMMON/PYSUBS/MSEL,MSUB(200),KFIN(2,-40:40),CKIN(200)
      COMMON/PYINT1/MINT(400),VINT(400)
      COMMON/PYINT2/ISET(200),KFPR(200,2),COEF(200,20),ICOL(40,4,2)
      COMMON/PYINT5/NGEN(0:200,3),XSEC(0:200,3)
      COMMON/PYINT9/DXSEC(0:200)
      DOUBLE PRECISION DXSEC
      SAVE /LUJETS/,/LUDAT1/
      SAVE /PYSUBS/,/PYPARS/,/PYINT1/,/PYINT2/,/PYINT5/,/PYINT9/
 
C...Calculate Monte Carlo estimates of cross-sections.
      ISUB=MINT(1)
      IF(MSTP(111).NE.-1) NGEN(ISUB,3)=NGEN(ISUB,3)+1
      NGEN(0,3)=NGEN(0,3)+1
      XSEC(0,3)=0.
      DO 100 I=1,200
      IF(I.EQ.96.OR.I.EQ.97) THEN
        XSEC(I,3)=0.
      ELSEIF(MSUB(95).EQ.1.AND.(I.EQ.11.OR.I.EQ.12.OR.I.EQ.13.OR.
     &I.EQ.28.OR.I.EQ.53.OR.I.EQ.68)) THEN
        XSEC(I,3)=DXSEC(96)*NGEN(I,3)/MAX(1.,FLOAT(NGEN(96,1))*
     &  FLOAT(NGEN(96,2)))
      ELSEIF(MSUB(I).EQ.0.OR.NGEN(I,1).EQ.0) THEN
        XSEC(I,3)=0.
      ELSEIF(NGEN(I,2).EQ.0) THEN
        XSEC(I,3)=DXSEC(I)*NGEN(0,3)/(FLOAT(NGEN(I,1))*
     &  FLOAT(NGEN(0,2)))
      ELSE
        XSEC(I,3)=DXSEC(I)*NGEN(I,3)/(FLOAT(NGEN(I,1))*
     &  FLOAT(NGEN(I,2)))
      ENDIF
      XSEC(0,3)=XSEC(0,3)+XSEC(I,3)
  100 CONTINUE
 
C...Rescale to known low-pT cross-section for standard QCD processes.
      IF(MSUB(95).EQ.1) THEN
        XSECH=XSEC(11,3)+XSEC(12,3)+XSEC(13,3)+XSEC(28,3)+XSEC(53,3)+
     &  XSEC(68,3)+XSEC(95,3)
        XSECW=DXSEC(97)/MAX(1.,FLOAT(NGEN(97,1)))
        IF(XSECH.GT.1E-10.AND.XSECW.GT.1E-10) THEN
          FAC=XSECW/XSECH
          XSEC(11,3)=FAC*XSEC(11,3)
          XSEC(12,3)=FAC*XSEC(12,3)
          XSEC(13,3)=FAC*XSEC(13,3)
          XSEC(28,3)=FAC*XSEC(28,3)
          XSEC(53,3)=FAC*XSEC(53,3)
          XSEC(68,3)=FAC*XSEC(68,3)
          XSEC(95,3)=FAC*XSEC(95,3)
          XSEC(0,3)=XSEC(0,3)-XSECH+XSECW
        ENDIF
      ENDIF
 
C...Save information for gamma-p and gamma-gamma.
      IF(MINT(121).GT.1) THEN
        IGA=MINT(122)
        CALL PYSAVE(2,IGA)
        CALL PYSAVE(5,0)
      ENDIF
 
C...Reset information on hard interaction.
      DO 110 J=1,200
      MSTI(J)=0
      PARI(J)=0.
  110 CONTINUE
 
C...Copy integer valued information from MINT into MSTI.
      DO 120 J=1,31
      MSTI(J)=MINT(J)
  120 CONTINUE
      IF(MINT(121).GT.1) MSTI(9)=MINT(122)
 
C...Store cross-section variables in PARI.
      PARI(1)=XSEC(0,3)
      PARI(2)=XSEC(0,3)/MINT(5)
      PARI(9)=VINT(99)
      PARI(10)=VINT(100)
      VINT(98)=VINT(98)+VINT(100)
      IF(MSTP(142).EQ.1) PARI(2)=XSEC(0,3)/VINT(98)
 
C...Store kinematics variables in PARI.
      PARI(11)=VINT(1)
      PARI(12)=VINT(2)
      IF(ISUB.NE.95) THEN
        DO 130 J=13,26
        PARI(J)=VINT(30+J)
  130   CONTINUE
        PARI(31)=VINT(141)
        PARI(32)=VINT(142)
        PARI(33)=VINT(41)
        PARI(34)=VINT(42)
        PARI(35)=PARI(33)-PARI(34)
        PARI(36)=VINT(21)
        PARI(37)=VINT(22)
        PARI(38)=VINT(26)
        PARI(39)=VINT(157)
        PARI(40)=VINT(158)
        PARI(41)=VINT(23)
        PARI(42)=2.*VINT(47)/VINT(1)
      ENDIF
 
C...Store information on scattered partons in PARI.
      IF(ISUB.NE.95.AND.MINT(7)*MINT(8).NE.0) THEN
        DO 140 IS=7,8
        I=MINT(IS)
        PARI(36+IS)=P(I,3)/VINT(1)
        PARI(38+IS)=P(I,4)/VINT(1)
        PR=MAX(1E-20,P(I,5)**2+P(I,1)**2+P(I,2)**2)
        PARI(40+IS)=SIGN(LOG(MIN((SQRT(PR+P(I,3)**2)+ABS(P(I,3)))/
     &  SQRT(PR),1E20)),P(I,3))
        PR=MAX(1E-20,P(I,1)**2+P(I,2)**2)
        PARI(42+IS)=SIGN(LOG(MIN((SQRT(PR+P(I,3)**2)+ABS(P(I,3)))/
     &  SQRT(PR),1E20)),P(I,3))
        PARI(44+IS)=P(I,3)/SQRT(1E-20+P(I,1)**2+P(I,2)**2+P(I,3)**2)
        PARI(46+IS)=ULANGL(P(I,3),SQRT(P(I,1)**2+P(I,2)**2))
        PARI(48+IS)=ULANGL(P(I,1),P(I,2))
  140   CONTINUE
      ENDIF
 
C...Store sum up transverse and longitudinal momenta.
      PARI(65)=2.*PARI(17)
      IF(ISUB.LE.90.OR.ISUB.GE.95) THEN
        DO 150 I=MSTP(126)+1,N
        IF(K(I,1).LE.0.OR.K(I,1).GT.10) GOTO 150
        PT=SQRT(P(I,1)**2+P(I,2)**2)
        PARI(69)=PARI(69)+PT
        IF(I.LE.MINT(52)) PARI(66)=PARI(66)+PT
        IF(I.GT.MINT(52).AND.I.LE.MINT(53)) PARI(68)=PARI(68)+PT
  150   CONTINUE
        PARI(67)=PARI(68)
        PARI(71)=VINT(151)
        PARI(72)=VINT(152)
        PARI(73)=VINT(151)
        PARI(74)=VINT(152)
      ELSE
        PARI(66)=PARI(65)
        PARI(69)=PARI(65)
      ENDIF
 
C...Store various other pieces of information into PARI.
      PARI(61)=VINT(148)
      PARI(75)=VINT(155)
      PARI(76)=VINT(156)
      PARI(77)=VINT(159)
      PARI(78)=VINT(160)
      PARI(81)=VINT(138)
 
C...Set information for LUTABU.
      IF(ISET(ISUB).EQ.1.OR.ISET(ISUB).EQ.3) THEN
        MSTU(161)=MINT(21)
        MSTU(162)=0
      ELSEIF(ISET(ISUB).EQ.5) THEN
        MSTU(161)=MINT(23)
        MSTU(162)=0
      ELSE
        MSTU(161)=MINT(21)
        MSTU(162)=MINT(22)
      ENDIF
 
      RETURN
      END
