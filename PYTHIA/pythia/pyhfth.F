 
C*********************************************************************
 
      FUNCTION PYHFTH(SH,SQM,FRATT)
 
C...Gives threshold attractive/repulsive factor for heavy flavour
C...production.
      COMMON/LUDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
      COMMON/PYPARS/MSTP(200),PARP(200),MSTI(200),PARI(200)
      COMMON/PYINT1/MINT(400),VINT(400)
      SAVE /LUDAT1/
      SAVE /PYPARS/,/PYINT1/
 
C...Value for alpha_strong.
      IF(MSTP(35).LE.1) THEN
        ALSSG=PARP(35)
      ELSE
        MST115=MSTU(115)
        MSTU(115)=MSTP(36)
        Q2BN=SQRT(MAX(1.,SQM*((SQRT(SH)-2.*SQRT(SQM))**2+PARP(36)**2)))
        ALSSG=ULALPS(Q2BN)
        MSTU(115)=MST115
      ENDIF
 
C...Evaluate attractive and repulsive factors.
      XATTR=4.*PARU(1)*ALSSG/(3.*SQRT(MAX(1E-20,1.-4.*SQM/SH)))
      FATTR=XATTR/(1.-EXP(-MIN(50.,XATTR)))
      XREPU=PARU(1)*ALSSG/(6.*SQRT(MAX(1E-20,1.-4.*SQM/SH)))
      FREPU=XREPU/(EXP(MIN(50.,XREPU))-1.)
      PYHFTH=FRATT*FATTR+(1.-FRATT)*FREPU
      VINT(138)=PYHFTH
 
      RETURN
      END
