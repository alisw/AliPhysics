*CMZ :          17/07/98  15.45.01  by  Federico Carminati
*-- Author :
C*********************************************************************

      SUBROUTINE LUEDIT(MEDIT)

C...Purpose: to perform global manipulations on the event record,
C...in particular to exclude unstable or undetectable partons/particles.
*KEEP,LUJETS.
      COMMON /LUJETS/ N,K(200000,5),P(200000,5),V(200000,5)
      SAVE /LUJETS/
*KEEP,LUDAT1.
      COMMON /LUDAT1/ MSTU(200),PARU(200),MSTJ(200),PARJ(200)
      SAVE /LUDAT1/
*KEEP,LUDAT2.
      COMMON /LUDAT2/ KCHG(500,3),PMAS(500,4),PARF(2000),VCKM(4,4)
      SAVE /LUDAT2/
*KEEP,SHWATE.
      COMMON /SHWATE/ WEI(200000)

*KEND.
      DIMENSION NS(2),PTS(2),PLS(2)

C...Remove unwanted partons/particles.
      IF((MEDIT.GE.0.AND.MEDIT.LE.3).OR.MEDIT.EQ.5) THEN
        IMAX=N
        IF(MSTU(2).GT.0) IMAX=MSTU(2)
        I1=MAX(1,MSTU(1))-1
        DO 110 I=MAX(1,MSTU(1)),IMAX
        IF(K(I,1).EQ.0.OR.K(I,1).GT.20) GOTO 110
        IF(MEDIT.EQ.1) THEN
          IF(K(I,1).GT.10) GOTO 110
        ELSEIF(MEDIT.EQ.2) THEN
          IF(K(I,1).GT.10) GOTO 110
          KC=LUCOMP(K(I,2))
          IF(KC.EQ.0.OR.KC.EQ.12.OR.KC.EQ.14.OR.KC.EQ.16.OR.KC.EQ.18)
     &    GOTO 110
        ELSEIF(MEDIT.EQ.3) THEN
          IF(K(I,1).GT.10) GOTO 110
          KC=LUCOMP(K(I,2))
          IF(KC.EQ.0) GOTO 110
          IF(KCHG(KC,2).EQ.0.AND.LUCHGE(K(I,2)).EQ.0) GOTO 110
        ELSEIF(MEDIT.EQ.5) THEN
          IF(K(I,1).EQ.13.OR.K(I,1).EQ.14) GOTO 110
          KC=LUCOMP(K(I,2))
          IF(KC.EQ.0) GOTO 110
          IF(K(I,1).GE.11.AND.KCHG(KC,2).EQ.0) GOTO 110
        ENDIF

C...Pack remaining partons/particles. Origin no longer known.
        I1=I1+1
        DO 100 J=1,5
        K(I1,J)=K(I,J)
        P(I1,J)=P(I,J)
  100   V(I1,J)=V(I,J)
        K(I1,3)=0
        WEI(I1)=WEI(I)
  110   CONTINUE
        IF(I1.LT.N) MSTU(3)=0
        IF(I1.LT.N) MSTU(70)=0
        N=I1

C...Selective removal of class of entries. New position of retained.
      ELSEIF(MEDIT.GE.11.AND.MEDIT.LE.15) THEN
        I1=0
        DO 120 I=1,N
        K(I,3)=MOD(K(I,3),MSTU(5))
        IF(MEDIT.EQ.11.AND.K(I,1).LT.0) GOTO 120
        IF(MEDIT.EQ.12.AND.K(I,1).EQ.0) GOTO 120
        IF(MEDIT.EQ.13.AND.(K(I,1).EQ.11.OR.K(I,1).EQ.12.OR.
     &  K(I,1).EQ.15).AND.K(I,2).NE.94) GOTO 120
        IF(MEDIT.EQ.14.AND.(K(I,1).EQ.13.OR.K(I,1).EQ.14.OR.
     &  K(I,2).EQ.94)) GOTO 120
        IF(MEDIT.EQ.15.AND.K(I,1).GE.21) GOTO 120
        I1=I1+1
        K(I,3)=K(I,3)+MSTU(5)*I1
  120   CONTINUE

C...Find new event history information and replace old.
        DO 140 I=1,N
        IF(K(I,1).LE.0.OR.K(I,1).GT.20.OR.K(I,3)/MSTU(5).EQ.0) GOTO 140
        ID=I
  130   IM=MOD(K(ID,3),MSTU(5))
        IF(MEDIT.EQ.13.AND.IM.GT.0.AND.IM.LE.N) THEN
          IF((K(IM,1).EQ.11.OR.K(IM,1).EQ.12.OR.K(IM,1).EQ.15).AND.
     &    K(IM,2).NE.94) THEN
            ID=IM
            GOTO 130
          ENDIF
        ELSEIF(MEDIT.EQ.14.AND.IM.GT.0.AND.IM.LE.N) THEN
          IF(K(IM,1).EQ.13.OR.K(IM,1).EQ.14.OR.K(IM,2).EQ.94) THEN
            ID=IM
            GOTO 130
          ENDIF
        ENDIF
        K(I,3)=MSTU(5)*(K(I,3)/MSTU(5))
        IF(IM.NE.0) K(I,3)=K(I,3)+K(IM,3)/MSTU(5)
        IF(K(I,1).NE.3.AND.K(I,1).NE.13.AND.K(I,1).NE.14) THEN
          IF(K(I,4).GT.0.AND.K(I,4).LE.MSTU(4)) K(I,4)=
     &    K(K(I,4),3)/MSTU(5)
          IF(K(I,5).GT.0.AND.K(I,5).LE.MSTU(4)) K(I,5)=
     &    K(K(I,5),3)/MSTU(5)
        ELSE
          KCM=MOD(K(I,4)/MSTU(5),MSTU(5))
          IF(KCM.GT.0.AND.KCM.LE.MSTU(4)) KCM=K(KCM,3)/MSTU(5)
          KCD=MOD(K(I,4),MSTU(5))
          IF(KCD.GT.0.AND.KCD.LE.MSTU(4)) KCD=K(KCD,3)/MSTU(5)
          K(I,4)=MSTU(5)**2*(K(I,4)/MSTU(5)**2)+MSTU(5)*KCM+KCD
          KCM=MOD(K(I,5)/MSTU(5),MSTU(5))
          IF(KCM.GT.0.AND.KCM.LE.MSTU(4)) KCM=K(KCM,3)/MSTU(5)
          KCD=MOD(K(I,5),MSTU(5))
          IF(KCD.GT.0.AND.KCD.LE.MSTU(4)) KCD=K(KCD,3)/MSTU(5)
          K(I,5)=MSTU(5)**2*(K(I,5)/MSTU(5)**2)+MSTU(5)*KCM+KCD
        ENDIF
  140   CONTINUE

C...Pack remaining entries.
        I1=0
        DO 160 I=1,N
        IF(K(I,3)/MSTU(5).EQ.0) GOTO 160
        I1=I1+1
        DO 150 J=1,5
        K(I1,J)=K(I,J)
        P(I1,J)=P(I,J)
  150   V(I1,J)=V(I,J)
        K(I1,3)=MOD(K(I1,3),MSTU(5))
  160   CONTINUE
        IF(I1.LT.N) MSTU(3)=0
        IF(I1.LT.N) MSTU(70)=0
        N=I1

C...Save top entries at bottom of LUJETS commonblock.
      ELSEIF(MEDIT.EQ.21) THEN
        IF(2*N.GE.MSTU(4)) THEN
          CALL LUERRM(11,'(LUEDIT:) no more memory left in LUJETS')
          RETURN
        ENDIF
        DO 170 I=1,N
        DO 170 J=1,5
        K(MSTU(4)-I,J)=K(I,J)
        P(MSTU(4)-I,J)=P(I,J)
  170   V(MSTU(4)-I,J)=V(I,J)
        MSTU(32)=N

C...Restore bottom entries of commonblock LUJETS to top.
      ELSEIF(MEDIT.EQ.22) THEN
        DO 180 I=1,MSTU(32)
        DO 180 J=1,5
        K(I,J)=K(MSTU(4)-I,J)
        P(I,J)=P(MSTU(4)-I,J)
  180   V(I,J)=V(MSTU(4)-I,J)
        N=MSTU(32)

C...Mark primary entries at top of commonblock LUJETS as untreated.
      ELSEIF(MEDIT.EQ.23) THEN
        I1=0
        DO 190 I=1,N
        KH=K(I,3)
        IF(KH.GE.1) THEN
          IF(K(KH,1).GT.20) KH=0
        ENDIF
        IF(KH.NE.0) GOTO 200
        I1=I1+1
  190   IF(K(I,1).GT.10.AND.K(I,1).LE.20) K(I,1)=K(I,1)-10
  200   N=I1

C...Place largest axis along z axis and second largest in xy plane.
      ELSEIF(MEDIT.EQ.31.OR.MEDIT.EQ.32) THEN
        CALL LUDBRB(1,N+MSTU(3),0.,-ULANGL(P(MSTU(61),1),
     &  P(MSTU(61),2)),0D0,0D0,0D0)
        CALL LUDBRB(1,N+MSTU(3),-ULANGL(P(MSTU(61),3),
     &  P(MSTU(61),1)),0.,0D0,0D0,0D0)
        CALL LUDBRB(1,N+MSTU(3),0.,-ULANGL(P(MSTU(61)+1,1),
     &  P(MSTU(61)+1,2)),0D0,0D0,0D0)
        IF(MEDIT.EQ.31) RETURN

C...Rotate to put slim jet along +z axis.
        DO 210 IS=1,2
        NS(IS)=0
        PTS(IS)=0.
  210   PLS(IS)=0.
        DO 220 I=1,N
        IF(K(I,1).LE.0.OR.K(I,1).GT.10) GOTO 220
        IF(MSTU(41).GE.2) THEN
          KC=LUCOMP(K(I,2))
          IF(KC.EQ.0.OR.KC.EQ.12.OR.KC.EQ.14.OR.KC.EQ.16.OR.
     &    KC.EQ.18) GOTO 220
          IF(MSTU(41).GE.3.AND.KCHG(KC,2).EQ.0.AND.LUCHGE(K(I,2)).EQ.0)
     &    GOTO 220
        ENDIF
        IS=2.-SIGN(0.5,P(I,3))
        NS(IS)=NS(IS)+1
        PTS(IS)=PTS(IS)+SQRT(P(I,1)**2+P(I,2)**2)
  220   CONTINUE
        IF(NS(1)*PTS(2)**2.LT.NS(2)*PTS(1)**2)
     &  CALL LUDBRB(1,N+MSTU(3),PARU(1),0.,0D0,0D0,0D0)

C...Rotate to put second largest jet into -z,+x quadrant.
        DO 230 I=1,N
        IF(P(I,3).GE.0.) GOTO 230
        IF(K(I,1).LE.0.OR.K(I,1).GT.10) GOTO 230
        IF(MSTU(41).GE.2) THEN
          KC=LUCOMP(K(I,2))
          IF(KC.EQ.0.OR.KC.EQ.12.OR.KC.EQ.14.OR.KC.EQ.16.OR.
     &    KC.EQ.18) GOTO 230
          IF(MSTU(41).GE.3.AND.KCHG(KC,2).EQ.0.AND.LUCHGE(K(I,2)).EQ.0)
     &    GOTO 230
        ENDIF
        IS=2.-SIGN(0.5,P(I,1))
        PLS(IS)=PLS(IS)-P(I,3)
  230   CONTINUE
        IF(PLS(2).GT.PLS(1)) CALL LUDBRB(1,N+MSTU(3),0.,PARU(1),
     &  0D0,0D0,0D0)
      ENDIF

      RETURN
      END
