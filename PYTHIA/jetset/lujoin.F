 
C********************************************************************* 
 
      SUBROUTINE LUJOIN(NJOIN,IJOIN) 
 
C...Purpose: to connect a sequence of partons with colour flow indices, 
C...as required for subsequent shower evolution (or other operations). 
      COMMON/LUJETS/N,K(4000,5),P(4000,5),V(4000,5) 
      COMMON/LUDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200) 
      COMMON/LUDAT2/KCHG(500,3),PMAS(500,4),PARF(2000),VCKM(4,4) 
      SAVE /LUJETS/,/LUDAT1/,/LUDAT2/ 
      DIMENSION IJOIN(*) 
 
C...Check that partons are of right types to be connected. 
      IF(NJOIN.LT.2) GOTO 120 
      KQSUM=0 
      DO 100 IJN=1,NJOIN 
      I=IJOIN(IJN) 
      IF(I.LE.0.OR.I.GT.N) GOTO 120 
      IF(K(I,1).LT.1.OR.K(I,1).GT.3) GOTO 120 
      KC=LUCOMP(K(I,2)) 
      IF(KC.EQ.0) GOTO 120 
      KQ=KCHG(KC,2)*ISIGN(1,K(I,2)) 
      IF(KQ.EQ.0) GOTO 120 
      IF(IJN.NE.1.AND.IJN.NE.NJOIN.AND.KQ.NE.2) GOTO 120 
      IF(KQ.NE.2) KQSUM=KQSUM+KQ 
      IF(IJN.EQ.1) KQS=KQ 
  100 CONTINUE 
      IF(KQSUM.NE.0) GOTO 120 
 
C...Connect the partons sequentially (closing for gluon loop). 
      KCS=(9-KQS)/2 
      IF(KQS.EQ.2) KCS=INT(4.5+RLU(0)) 
      DO 110 IJN=1,NJOIN 
      I=IJOIN(IJN) 
      K(I,1)=3 
      IF(IJN.NE.1) IP=IJOIN(IJN-1) 
      IF(IJN.EQ.1) IP=IJOIN(NJOIN) 
      IF(IJN.NE.NJOIN) IN=IJOIN(IJN+1) 
      IF(IJN.EQ.NJOIN) IN=IJOIN(1) 
      K(I,KCS)=MSTU(5)*IN 
      K(I,9-KCS)=MSTU(5)*IP 
      IF(IJN.EQ.1.AND.KQS.NE.2) K(I,9-KCS)=0 
      IF(IJN.EQ.NJOIN.AND.KQS.NE.2) K(I,KCS)=0 
  110 CONTINUE 
 
C...Error exit: no action taken. 
      RETURN 
  120 CALL LUERRM(12, 
     &'(LUJOIN:) given entries can not be joined by one string') 
 
      RETURN 
      END 
