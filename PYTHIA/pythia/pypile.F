 
C*********************************************************************
 
      SUBROUTINE PYPILE(MPILE)
 
C...Initializes multiplicity distribution and selects mutliplicity
C...of pileup events, i.e. several events occuring at the same
C...beam crossing.
      COMMON/LUDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
      COMMON/PYPARS/MSTP(200),PARP(200),MSTI(200),PARI(200)
      COMMON/PYINT1/MINT(400),VINT(400)
      COMMON/PYINT7/SIGT(0:6,0:6,0:5)
      SAVE /LUDAT1/
      SAVE /PYPARS/,/PYINT1/,/PYINT7/
      DIMENSION WTI(0:200)
      SAVE IMIN,IMAX,WTI,WTS
 
C...Sum of allowed cross-sections for pileup events.
      IF(MPILE.EQ.1) THEN
        VINT(131)=SIGT(0,0,5)
        IF(MSTP(132).GE.2) VINT(131)=VINT(131)+SIGT(0,0,4)
        IF(MSTP(132).GE.3) VINT(131)=VINT(131)+SIGT(0,0,2)+SIGT(0,0,3)
        IF(MSTP(132).GE.4) VINT(131)=VINT(131)+SIGT(0,0,1)
        IF(MSTP(133).LE.0) RETURN
 
C...Initialize multiplicity distribution at maximum.
        XNAVE=VINT(131)*PARP(131)
        IF(XNAVE.GT.120.) WRITE(MSTU(11),5000) XNAVE
        INAVE=MAX(1,MIN(200,NINT(XNAVE)))
        WTI(INAVE)=1.
        WTS=WTI(INAVE)
        WTN=WTI(INAVE)*INAVE
 
C...Find shape of multiplicity distribution below maximum.
        IMIN=INAVE
        DO 100 I=INAVE-1,1,-1
        IF(MSTP(133).EQ.1) WTI(I)=WTI(I+1)*(I+1)/XNAVE
        IF(MSTP(133).GE.2) WTI(I)=WTI(I+1)*I/XNAVE
        IF(WTI(I).LT.1E-6) GOTO 110
        WTS=WTS+WTI(I)
        WTN=WTN+WTI(I)*I
        IMIN=I
  100   CONTINUE
 
C...Find shape of multiplicity distribution above maximum.
  110   IMAX=INAVE
        DO 120 I=INAVE+1,200
        IF(MSTP(133).EQ.1) WTI(I)=WTI(I-1)*XNAVE/I
        IF(MSTP(133).GE.2) WTI(I)=WTI(I-1)*XNAVE/(I-1)
        IF(WTI(I).LT.1E-6) GOTO 130
        WTS=WTS+WTI(I)
        WTN=WTN+WTI(I)*I
        IMAX=I
  120   CONTINUE
  130   VINT(132)=XNAVE
        VINT(133)=WTN/WTS
        IF(MSTP(133).EQ.1.AND.IMIN.EQ.1) VINT(134)=
     &  WTS/(WTS+WTI(1)/XNAVE)
        IF(MSTP(133).EQ.1.AND.IMIN.GT.1) VINT(134)=1.
        IF(MSTP(133).GE.2) VINT(134)=XNAVE
 
C...Pick multiplicity of pileup events.
      ELSE
        IF(MSTP(133).LE.0) THEN
          MINT(81)=MAX(1,MSTP(134))
        ELSE
          WTR=WTS*RLU(0)
          DO 140 I=IMIN,IMAX
          MINT(81)=I
          WTR=WTR-WTI(I)
          IF(WTR.LE.0.) GOTO 150
  140     CONTINUE
  150     CONTINUE
        ENDIF
      ENDIF
 
C...Format statement for error message.
 5000 FORMAT(1X,'Warning: requested average number of events per bunch',
     &'crossing too large, ',1P,E12.4)
 
      RETURN
      END
