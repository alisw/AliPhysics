*CMZ :          17/07/98  15.44.34  by  Federico Carminati
*-- Author :
C*********************************************************************

      FUNCTION PLU(I,J)

C...Purpose: to provide various real-valued event related data.
*KEEP,LUJETS.
      COMMON /LUJETS/ N,K(200000,5),P(200000,5),V(200000,5)
      SAVE /LUJETS/
*KEEP,LUDAT1.
      COMMON /LUDAT1/ MSTU(200),PARU(200),MSTJ(200),PARJ(200)
      SAVE /LUDAT1/
*KEEP,LUDAT2.
      COMMON /LUDAT2/ KCHG(500,3),PMAS(500,4),PARF(2000),VCKM(4,4)
      SAVE /LUDAT2/
*KEND.
      DIMENSION PSUM(4)

C...Set default value. For I = 0 sum of momenta or charges,
C...or invariant mass of system.
      PLU=0.
      IF(I.LT.0.OR.I.GT.MSTU(4).OR.J.LE.0) THEN
      ELSEIF(I.EQ.0.AND.J.LE.4) THEN
        DO 100 I1=1,N
  100   IF(K(I1,1).GT.0.AND.K(I1,1).LE.10) PLU=PLU+P(I1,J)
      ELSEIF(I.EQ.0.AND.J.EQ.5) THEN
        DO 110 J1=1,4
        PSUM(J1)=0.
        DO 110 I1=1,N
  110   IF(K(I1,1).GT.0.AND.K(I1,1).LE.10) PSUM(J1)=PSUM(J1)+P(I1,J1)
        PLU=SQRT(MAX(0.,PSUM(4)**2-PSUM(1)**2-PSUM(2)**2-PSUM(3)**2))
      ELSEIF(I.EQ.0.AND.J.EQ.6) THEN
        DO 120 I1=1,N
  120   IF(K(I1,1).GT.0.AND.K(I1,1).LE.10) PLU=PLU+LUCHGE(K(I1,2))/3.
      ELSEIF(I.EQ.0) THEN

C...Direct readout of P matrix.
      ELSEIF(J.LE.5) THEN
        PLU=P(I,J)

C...Charge, total momentum, transverse momentum, transverse mass.
      ELSEIF(J.LE.12) THEN
        IF(J.EQ.6) PLU=LUCHGE(K(I,2))/3.
        IF(J.EQ.7.OR.J.EQ.8) PLU=P(I,1)**2+P(I,2)**2+P(I,3)**2
        IF(J.EQ.9.OR.J.EQ.10) PLU=P(I,1)**2+P(I,2)**2
        IF(J.EQ.11.OR.J.EQ.12) PLU=P(I,5)**2+P(I,1)**2+P(I,2)**2
        IF(J.EQ.8.OR.J.EQ.10.OR.J.EQ.12) PLU=SQRT(PLU)

C...Theta and phi angle in radians or degrees.
      ELSEIF(J.LE.16) THEN
        IF(J.LE.14) PLU=ULANGL(P(I,3),SQRT(P(I,1)**2+P(I,2)**2))
        IF(J.GE.15) PLU=ULANGL(P(I,1),P(I,2))
        IF(J.EQ.14.OR.J.EQ.16) PLU=PLU*180./PARU(1)

C...True rapidity, rapidity with pion mass, pseudorapidity.
      ELSEIF(J.LE.19) THEN
        PMR=0.
        IF(J.EQ.17) PMR=P(I,5)
        IF(J.EQ.18) PMR=ULMASS(211)
        PR=MAX(1E-20,PMR**2+P(I,1)**2+P(I,2)**2)
        PLU=SIGN(LOG(MIN((SQRT(PR+P(I,3)**2)+ABS(P(I,3)))/SQRT(PR),
     &  1E20)),P(I,3))

C...Energy and momentum fractions (only to be used in CM frame).
      ELSEIF(J.LE.25) THEN
        IF(J.EQ.20) PLU=2.*SQRT(P(I,1)**2+P(I,2)**2+P(I,3)**2)/PARU(21)
        IF(J.EQ.21) PLU=2.*P(I,3)/PARU(21)
        IF(J.EQ.22) PLU=2.*SQRT(P(I,1)**2+P(I,2)**2)/PARU(21)
        IF(J.EQ.23) PLU=2.*P(I,4)/PARU(21)
        IF(J.EQ.24) PLU=(P(I,4)+P(I,3))/PARU(21)
        IF(J.EQ.25) PLU=(P(I,4)-P(I,3))/PARU(21)
      ENDIF

      RETURN
      END
