 
C***********************************************************************
 
      SUBROUTINE PYI3AU(EPS,RAT,Y3RE,Y3IM)
 
C...Calculates real and imaginary parts of the auxiliary function I3;
C...see R. K. Ellis, I. Hinchliffe, M. Soldate and J. J. van der Bij,
C...Nucl. Phys. B297 (1988) 221.
      COMMON/LUDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
      SAVE /LUDAT1/
 
      BE=0.5*(1.+SQRT(1.+RAT*EPS))
      IF(EPS.LT.1.) GA=0.5*(1.+SQRT(1.-EPS))
 
      IF(EPS.LT.0.) THEN
        IF(ABS(EPS).LT.1.E-4.AND.ABS(RAT*EPS).LT.1.E-4) THEN
          F3RE=PYSPEN(-0.25*EPS/(1.+0.25*(RAT-1.)*EPS),0.,1)-
     &    PYSPEN((1.-0.25*EPS)/(1.+0.25*(RAT-1.)*EPS),0.,1)+
     &    PYSPEN(0.25*(RAT+1.)*EPS/(1.+0.25*RAT*EPS),0.,1)-
     &    PYSPEN((RAT+1.)/RAT,0.,1)+0.5*(LOG(1.+0.25*RAT*EPS)**2-
     &    LOG(0.25*RAT*EPS)**2)+LOG(1.-0.25*EPS)*
     &    LOG((1.+0.25*(RAT-1.)*EPS)/(1.+0.25*RAT*EPS))+
     &    LOG(-0.25*EPS)*LOG(0.25*RAT*EPS/(1.+0.25*(RAT-1.)*EPS))
        ELSEIF(ABS(EPS).LT.1.E-4.AND.ABS(RAT*EPS).GE.1.E-4) THEN
          F3RE=PYSPEN(-0.25*EPS/(BE-0.25*EPS),0.,1)-
     &    PYSPEN((1.-0.25*EPS)/(BE-0.25*EPS),0.,1)+
     &    PYSPEN((BE-1.+0.25*EPS)/BE,0.,1)-
     &    PYSPEN((BE-1.+0.25*EPS)/(BE-1.),0.,1)+
     &    0.5*(LOG(BE)**2-LOG(BE-1.)**2)+
     &    LOG(1.-0.25*EPS)*LOG((BE-0.25*EPS)/BE)+
     &    LOG(-0.25*EPS)*LOG((BE-1.)/(BE-0.25*EPS))
        ELSEIF(ABS(EPS).GE.1.E-4.AND.ABS(RAT*EPS).LT.1.E-4) THEN
          F3RE=PYSPEN((GA-1.)/(GA+0.25*RAT*EPS),0.,1)-
     &    PYSPEN(GA/(GA+0.25*RAT*EPS),0.,1)+
     &    PYSPEN((1.+0.25*RAT*EPS-GA)/(1.+0.25*RAT*EPS),0.,1)-
     &    PYSPEN((1.+0.25*RAT*EPS-GA)/(0.25*RAT*EPS),0.,1)+
     &    0.5*(LOG(1.+0.25*RAT*EPS)**2-LOG(0.25*RAT*EPS)**2)+
     &    LOG(GA)*LOG((GA+0.25*RAT*EPS)/(1.+0.25*RAT*EPS))+
     &    LOG(GA-1.)*LOG(0.25*RAT*EPS/(GA+0.25*RAT*EPS))
        ELSE
          F3RE=PYSPEN((GA-1.)/(GA+BE-1.),0.,1)-
     &    PYSPEN(GA/(GA+BE-1.),0.,1)+PYSPEN((BE-GA)/BE,0.,1)-
     &    PYSPEN((BE-GA)/(BE-1.),0.,1)+0.5*(LOG(BE)**2-LOG(BE-1.)**2)+
     &    LOG(GA)*LOG((GA+BE-1.)/BE)+LOG(GA-1.)*LOG((BE-1.)/(GA+BE-1.))
        ENDIF
        F3IM=0.
      ELSEIF(EPS.LT.1.) THEN
        IF(ABS(EPS).LT.1.E-4.AND.ABS(RAT*EPS).LT.1.E-4) THEN
          F3RE=PYSPEN(-0.25*EPS/(1.+0.25*(RAT-1.)*EPS),0.,1)-
     &    PYSPEN((1.-0.25*EPS)/(1.+0.25*(RAT-1.)*EPS),0.,1)+
     &    PYSPEN((1.-0.25*EPS)/(-0.25*(RAT+1.)*EPS),0.,1)-
     &    PYSPEN(1./(RAT+1.),0.,1)+LOG((1.-0.25*EPS)/(0.25*EPS))*
     &    LOG((1.+0.25*(RAT-1.)*EPS)/(0.25*(RAT+1.)*EPS))
          F3IM=-PARU(1)*LOG((1.+0.25*(RAT-1.)*EPS)/(0.25*(RAT+1.)*EPS))
        ELSEIF(ABS(EPS).LT.1.E-4.AND.ABS(RAT*EPS).GE.1.E-4) THEN
          F3RE=PYSPEN(-0.25*EPS/(BE-0.25*EPS),0.,1)-
     &    PYSPEN((1.-0.25*EPS)/(BE-0.25*EPS),0.,1)+
     &    PYSPEN((1.-0.25*EPS)/(1.-0.25*EPS-BE),0.,1)-
     &    PYSPEN(-0.25*EPS/(1.-0.25*EPS-BE),0.,1)+
     &    LOG((1.-0.25*EPS)/(0.25*EPS))*
     &    LOG((BE-0.25*EPS)/(BE-1.+0.25*EPS))
          F3IM=-PARU(1)*LOG((BE-0.25*EPS)/(BE-1.+0.25*EPS))
        ELSEIF(ABS(EPS).GE.1.E-4.AND.ABS(RAT*EPS).LT.1.E-4) THEN
          F3RE=PYSPEN((GA-1.)/(GA+0.25*RAT*EPS),0.,1)-
     &    PYSPEN(GA/(GA+0.25*RAT*EPS),0.,1)+
     &    PYSPEN(GA/(GA-1.-0.25*RAT*EPS),0.,1)-
     &    PYSPEN((GA-1.)/(GA-1.-0.25*RAT*EPS),0.,1)+
     &    LOG(GA/(1.-GA))*LOG((GA+0.25*RAT*EPS)/(1.+0.25*RAT*EPS-GA))
          F3IM=-PARU(1)*LOG((GA+0.25*RAT*EPS)/(1.+0.25*RAT*EPS-GA))
        ELSE
          F3RE=PYSPEN((GA-1.)/(GA+BE-1.),0.,1)-
     &    PYSPEN(GA/(GA+BE-1.),0.,1)+PYSPEN(GA/(GA-BE),0.,1)-
     &    PYSPEN((GA-1.)/(GA-BE),0.,1)+LOG(GA/(1.-GA))*
     &    LOG((GA+BE-1.)/(BE-GA))
          F3IM=-PARU(1)*LOG((GA+BE-1.)/(BE-GA))
         ENDIF
      ELSE
        RSQ=EPS/(EPS-1.+(2.*BE-1.)**2)
        RCTHE=RSQ*(1.-2.*BE/EPS)
        RSTHE=SQRT(MAX(0.,RSQ-RCTHE**2))
        RCPHI=RSQ*(1.+2.*(BE-1.)/EPS)
        RSPHI=SQRT(MAX(0.,RSQ-RCPHI**2))
        R=SQRT(RSQ)
        THE=ACOS(MAX(-0.999999,MIN(0.999999,RCTHE/R)))
        PHI=ACOS(MAX(-0.999999,MIN(0.999999,RCPHI/R)))
        F3RE=PYSPEN(RCTHE,RSTHE,1)+PYSPEN(RCTHE,-RSTHE,1)-
     &  PYSPEN(RCPHI,RSPHI,1)-PYSPEN(RCPHI,-RSPHI,1)+
     &  (PHI-THE)*(PHI+THE-PARU(1))
        F3IM=PYSPEN(RCTHE,RSTHE,2)+PYSPEN(RCTHE,-RSTHE,2)-
     &  PYSPEN(RCPHI,RSPHI,2)-PYSPEN(RCPHI,-RSPHI,2)
      ENDIF
 
      Y3RE=2./(2.*BE-1.)*F3RE
      Y3IM=2./(2.*BE-1.)*F3IM
 
      RETURN
      END
