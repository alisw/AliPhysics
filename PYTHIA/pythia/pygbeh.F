 
C*********************************************************************
 
      SUBROUTINE PYGBEH(KF,X,Q2,P2,PM2,XPBH)
C...Purpose: to evaluate the Bethe-Heitler cross section for
C...heavy flavour production.
      DATA AEM2PI/0.0011614/
 
C...Reset output.
      XPBH=0.
      SIGBH=0.
 
C...Check kinematics limits.
      IF(X.GE.Q2/(4.*PM2+Q2+P2)) RETURN
      W2=Q2*(1.-X)/X-P2
      BETA2=1.-4.*PM2/W2
      IF(BETA2.LT.1E-10) RETURN
      RMQ=4.*PM2/Q2
 
C...Simple case: P2 = 0.
      IF(P2.LT.1E-4) THEN
        BETA=SQRT(BETA2)
        IF(BETA.LT.0.99) THEN
          XBL=LOG((1.+BETA)/(1.-BETA))
        ELSE
          XBL=LOG((1.+BETA)**2*W2/(4.*PM2))
        ENDIF
        SIGBH=BETA*(8.*X*(1.-X)-1.-RMQ*X*(1.-X))+
     &  XBL*(X**2+(1.-X)**2+RMQ*X*(1.-3.*X)-0.5*RMQ**2*X**2)
 
C...Complicated case: P2 > 0, based on approximation of
C...C.T. Hill and G.G. Ross, Nucl. Phys. B148 (1979) 373
      ELSE
        RPQ=1.-4.*X**2*P2/Q2
        IF(RPQ.GT.1E-10) THEN
          RPBE=SQRT(RPQ*BETA2)
          IF(RPBE.LT.0.99) THEN
            XBL=LOG((1.+RPBE)/(1.-RPBE))
            XBI=2.*RPBE/(1.-RPBE**2)
          ELSE
            RPBESN=4.*PM2/W2+(4.*X**2*P2/Q2)*BETA2
            XBL=LOG((1.+RPBE)**2/RPBESN)
            XBI=2.*RPBE/RPBESN
          ENDIF
          SIGBH=BETA*(6.*X*(1.-X)-1.)+
     &    XBL*(X**2+(1.-X)**2+RMQ*X*(1.-3.*X)-0.5*RMQ**2*X**2)+
     &    XBI*(2.*X/Q2)*(PM2*X*(2.-RMQ)-P2*X)
        ENDIF
      ENDIF
 
C...Multiply by charge-squared etc. to get parton distribution.
      CHSQ=1./9.
      IF(IABS(KF).EQ.2.OR.IABS(KF).EQ.4) CHSQ=4./9.
      XPBH=3.*CHSQ*AEM2PI*X*SIGBH
 
      RETURN
      END
