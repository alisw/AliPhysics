 
C*********************************************************************
 
      SUBROUTINE PYFRAM(IFRAME)
 
C...Performs transformations between different coordinate frames.
      COMMON/LUDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
      COMMON/PYPARS/MSTP(200),PARP(200),MSTI(200),PARI(200)
      COMMON/PYINT1/MINT(400),VINT(400)
      SAVE /LUDAT1/
      SAVE /PYPARS/,/PYINT1/
 
C...Check that transformation can and should be done.
      IF(IFRAME.EQ.1.OR.IFRAME.EQ.2.OR.(IFRAME.EQ.3.AND.
     &MINT(91).EQ.1)) THEN
        IF(IFRAME.EQ.MINT(6)) RETURN
      ELSE
        WRITE(MSTU(11),5000) IFRAME,MINT(6)
        RETURN
      ENDIF
 
      IF(MINT(6).EQ.1) THEN
C...Transform from fixed target or user specified frame to
C...overall CM frame.
        CALL LUROBO(0.,0.,-VINT(8),-VINT(9),-VINT(10))
        CALL LUROBO(0.,-VINT(7),0.,0.,0.)
        CALL LUROBO(-VINT(6),0.,0.,0.,0.)
      ELSEIF(MINT(6).EQ.3) THEN
C...Transform from hadronic CM frame in DIS to overall CM frame.
        CALL LUROBO(-VINT(221),-VINT(222),-VINT(223),-VINT(224),
     &  -VINT(225))
      ENDIF
 
      IF(IFRAME.EQ.1) THEN
C...Transform from overall CM frame to fixed target or user specified
C...frame.
        CALL LUROBO(VINT(6),VINT(7),VINT(8),VINT(9),VINT(10))
      ELSEIF(IFRAME.EQ.3) THEN
C...Transform from overall CM frame to hadronic CM frame in DIS.
        CALL LUROBO(0.,0.,VINT(223),VINT(224),VINT(225))
        CALL LUROBO(0.,VINT(222),0.,0.,0.)
        CALL LUROBO(VINT(221),0.,0.,0.,0.)
      ENDIF
 
C...Set information about new frame.
      MINT(6)=IFRAME
      MSTI(6)=IFRAME
 
 5000 FORMAT(1X,'Error: illegal values in subroutine PYFRAM.',1X,
     &'No transformation performed.'/1X,'IFRAME =',1X,I5,'; MINT(6) =',
     &1X,I5)
 
      RETURN
      END
