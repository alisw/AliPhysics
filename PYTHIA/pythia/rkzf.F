 
*==================================================================
 
      FUNCTION RKZF(L1,P1,Q,RMB,L2,P2)
* COMPUTES THE SCALAR STRUCTURE
* U_BAR(L1,P1)(SLASH(Q)+RMB)U(L2,P2)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      COMPLEX RKZF,RKZPR,RKZSF
      COMPLEX ANSF(-1:1,1:4,1:8,-1:1,1:4)
      INTEGER DONF(-1:1,1:4,1:8,-1:1,1:4)
      COMMON/RKZFCO/ANSF,DONF
      DIMENSION P1(0:4),P2(0:4),Q(0:4),R(0:4)
* CHECK ON CORRECT LABEL INPUT
      IP1=IDINT(P1(4))
      IQ=IDINT(Q(4))
      IP2=IDINT(P2(4))
      IF(IABS(L1).NE.1.OR.IABS(L2).NE.1.OR.
     . IP1.LT.1.OR.IP1.GT.4            .OR.
     . IQ.LT.1.OR.IQ.GT.8              .OR.
     . IP2.LT.1.OR.IP2.GT.4) THEN
        WRITE(6,*) ' RKZF LABEL ERROR'
        WRITE(6,*) 'L1=',L1,' IP1=',IP1,' IQ=',IQ,
     .             ' L2=',L2,' IP2=',IP2
        STOP
      ENDIF
* CHECK WHETHER THIS ONE HAS BEEN CALCULATED ALREADY
      IF(DONF(L1,IP1,IQ,L2,IP2).EQ.0) THEN
* THIS ONE NOT DONE YET: DO IT AND STORE THE RESULT IN ARRAY 'ANSF'
        IF(L1.EQ.L2) THEN
          A=2D0*RKDOT(Q,P2)
C         IF(DABS(A).LT.(1D-10*P2(0)*Q(0))) THEN
C...The check above is extended to following.
          IF(ABS(A).LT.MAX(1D-8,ABS(1D-10*P2(0)*Q(0)))) THEN
            ANSF(L1,IP1,IQ,L2,IP2)=(0.,0.)
          ELSE
            A=RKDOT(Q,Q)/A
            DO 100 K=0,3
              R(K)=Q(K)-A*P2(K)
  100       CONTINUE
            IF(R(0).GT.0D0) THEN
              C=1D0
            ELSE
              DO 110 K=0,3
                R(K)=-R(K)
  110         CONTINUE
              C=-1D0
            ENDIF
            ANSF(L1,IP1,IQ,L2,IP2)=C*RKZPR(L1,P1,R)*RKZPR(-L1,R,P2)
          ENDIF
        ELSEIF(L1.EQ.-L2) THEN
          ANSF(L1,IP1,IQ,L2,IP2)=RMB*RKZSF(L1,P1,L2,P2)
        ELSE
          WRITE(6,*) ' ERROR IN RKZF: L1=',L1,'  L2=',L2
          STOP
        ENDIF
        RKZF=ANSF(L1,IP1,IQ,L2,IP2)
        DONF(L1,IP1,IQ,L2,IP2)=1
      ELSE
        RKZF=ANSF(L1,IP1,IQ,L2,IP2)
      ENDIF
      END
