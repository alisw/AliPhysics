 
C********************************************************************* 
 
      FUNCTION ULALPS(Q2) 
 
C...Purpose: to give the value of alpha_strong. 
      COMMON/LUDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200) 
      COMMON/LUDAT2/KCHG(500,3),PMAS(500,4),PARF(2000),VCKM(4,4) 
      SAVE /LUDAT1/,/LUDAT2/ 
 
C...Constant alpha_strong trivial. 
      IF(MSTU(111).LE.0) THEN 
        ULALPS=PARU(111) 
        MSTU(118)=MSTU(112) 
        PARU(117)=0. 
        PARU(118)=PARU(111) 
        RETURN 
      ENDIF 
 
C...Find effective Q2, number of flavours and Lambda. 
      Q2EFF=Q2 
      IF(MSTU(115).GE.2) Q2EFF=MAX(Q2,PARU(114)) 
      NF=MSTU(112) 
      ALAM2=PARU(112)**2 
  100 IF(NF.GT.MAX(2,MSTU(113))) THEN 
        Q2THR=PARU(113)*PMAS(NF,1)**2 
        IF(Q2EFF.LT.Q2THR) THEN 
          NF=NF-1 
          ALAM2=ALAM2*(Q2THR/ALAM2)**(2./(33.-2.*NF)) 
          GOTO 100 
        ENDIF 
      ENDIF 
  110 IF(NF.LT.MIN(8,MSTU(114))) THEN 
        Q2THR=PARU(113)*PMAS(NF+1,1)**2 
        IF(Q2EFF.GT.Q2THR) THEN 
          NF=NF+1 
          ALAM2=ALAM2*(ALAM2/Q2THR)**(2./(33.-2.*NF)) 
          GOTO 110 
        ENDIF 
      ENDIF 
      IF(MSTU(115).EQ.1) Q2EFF=Q2EFF+ALAM2 
      PARU(117)=SQRT(ALAM2) 
 
C...Evaluate first or second order alpha_strong. 
      B0=(33.-2.*NF)/6. 
      ALGQ=LOG(MAX(1.0001,Q2EFF/ALAM2)) 
      IF(MSTU(111).EQ.1) THEN 
        ULALPS=MIN(PARU(115),PARU(2)/(B0*ALGQ)) 
      ELSE 
        B1=(153.-19.*NF)/6. 
        ULALPS=MIN(PARU(115),PARU(2)/(B0*ALGQ)*(1.-B1*LOG(ALGQ)/ 
     &  (B0**2*ALGQ))) 
      ENDIF 
      MSTU(118)=NF 
      PARU(118)=ULALPS 
 
      RETURN 
      END 
