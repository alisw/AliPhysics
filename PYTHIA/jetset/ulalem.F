 
C********************************************************************* 
 
      FUNCTION ULALEM(Q2) 
 
C...Purpose: to calculate the running alpha_electromagnetic. 
      COMMON/LUDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200) 
      SAVE /LUDAT1/ 
 
C...Calculate real part of photon vacuum polarization. 
C...For leptons simplify by using asymptotic (Q^2 >> m^2) expressions. 
C...For hadrons use parametrization of H. Burkhardt et al. 
C...See R. Kleiss et al, CERN 89-08, vol. 3, pp. 129-131. 
      AEMPI=PARU(101)/(3.*PARU(1)) 
      IF(MSTU(101).LE.0.OR.Q2.LT.2E-6) THEN 
        RPIGG=0. 
      ELSEIF(MSTU(101).EQ.2.AND.Q2.LT.PARU(104)) THEN
        RPIGG=0.
      ELSEIF(MSTU(101).EQ.2) THEN
        RPIGG=1.-PARU(101)/PARU(103) 
      ELSEIF(Q2.LT.0.09) THEN 
        RPIGG=AEMPI*(13.4916+LOG(Q2))+0.00835*LOG(1.+Q2) 
      ELSEIF(Q2.LT.9.) THEN 
        RPIGG=AEMPI*(16.3200+2.*LOG(Q2))+0.00238*LOG(1.+3.927*Q2) 
      ELSEIF(Q2.LT.1E4) THEN 
        RPIGG=AEMPI*(13.4955+3.*LOG(Q2))+0.00165+0.00299*LOG(1.+Q2) 
      ELSE 
        RPIGG=AEMPI*(13.4955+3.*LOG(Q2))+0.00221+0.00293*LOG(1.+Q2) 
      ENDIF 
 
C...Calculate running alpha_em. 
      ULALEM=PARU(101)/(1.-RPIGG) 
      PARU(108)=ULALEM 
 
      RETURN 
      END 
