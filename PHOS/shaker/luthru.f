*CMZ :          17/07/98  15.44.34  by  Federico Carminati
*-- Author :
C*********************************************************************

      SUBROUTINE LUTHRU(THR,OBL)

C...Purpose: to perform thrust analysis to give thrust, oblateness
C...and the related event axes.
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
      DIMENSION TDI(3),TPR(3)

C...Take copy of particles that are to be considered in thrust analysis.
      NP=0
      PS=0.
      DO 100 I=1,N
      IF(K(I,1).LE.0.OR.K(I,1).GT.10) GOTO 100
      IF(MSTU(41).GE.2) THEN
        KC=LUCOMP(K(I,2))
        IF(KC.EQ.0.OR.KC.EQ.12.OR.KC.EQ.14.OR.KC.EQ.16.OR.
     &  KC.EQ.18) GOTO 100
        IF(MSTU(41).GE.3.AND.KCHG(KC,2).EQ.0.AND.LUCHGE(K(I,2)).EQ.0)
     &  GOTO 100
      ENDIF
      IF(N+NP+MSTU(44)+15.GE.MSTU(4)-MSTU(32)-5) THEN
        CALL LUERRM(11,'(LUTHRU:) no more memory left in LUJETS')
        THR=-2.
        OBL=-2.
        RETURN
      ENDIF
      NP=NP+1
      K(N+NP,1)=23
      P(N+NP,1)=P(I,1)
      P(N+NP,2)=P(I,2)
      P(N+NP,3)=P(I,3)
      P(N+NP,4)=SQRT(P(I,1)**2+P(I,2)**2+P(I,3)**2)
      P(N+NP,5)=1.
      IF(ABS(PARU(42)-1.).GT.0.001) P(N+NP,5)=P(N+NP,4)**(PARU(42)-1.)
      PS=PS+P(N+NP,4)*P(N+NP,5)
  100 CONTINUE

C...Very low multiplicities (0 or 1) not considered.
      IF(NP.LE.1) THEN
        CALL LUERRM(8,'(LUTHRU:) too few particles for analysis')
        THR=-1.
        OBL=-1.
        RETURN
      ENDIF

C...Loop over thrust and major. T axis along z direction in latter case.
      DO 280 ILD=1,2
      IF(ILD.EQ.2) THEN
        K(N+NP+1,1)=31
        PHI=ULANGL(P(N+NP+1,1),P(N+NP+1,2))
        MSTU(33)=1
        CALL LUDBRB(N+1,N+NP+1,0.,-PHI,0D0,0D0,0D0)
        THE=ULANGL(P(N+NP+1,3),P(N+NP+1,1))
        CALL LUDBRB(N+1,N+NP+1,-THE,0.,0D0,0D0,0D0)
      ENDIF

C...Find and order particles with highest p (pT for major).
      DO 110 ILF=N+NP+4,N+NP+MSTU(44)+4
  110 P(ILF,4)=0.
      DO 150 I=N+1,N+NP
      IF(ILD.EQ.2) P(I,4)=SQRT(P(I,1)**2+P(I,2)**2)
      DO 120 ILF=N+NP+MSTU(44)+3,N+NP+4,-1
      IF(P(I,4).LE.P(ILF,4)) GOTO 130
      DO 120 J=1,5
  120 P(ILF+1,J)=P(ILF,J)
      ILF=N+NP+3
  130 DO 140 J=1,5
  140 P(ILF+1,J)=P(I,J)
  150 CONTINUE

C...Find and order initial axes with highest thrust (major).
      DO 160 ILG=N+NP+MSTU(44)+5,N+NP+MSTU(44)+15
  160 P(ILG,4)=0.
      NC=2**(MIN(MSTU(44),NP)-1)
      DO 220 ILC=1,NC
      DO 170 J=1,3
  170 TDI(J)=0.
      DO 180 ILF=1,MIN(MSTU(44),NP)
      SGN=P(N+NP+ILF+3,5)
      IF(2**ILF*((ILC+2**(ILF-1)-1)/2**ILF).GE.ILC) SGN=-SGN
      DO 180 J=1,4-ILD
  180 TDI(J)=TDI(J)+SGN*P(N+NP+ILF+3,J)
      TDS=TDI(1)**2+TDI(2)**2+TDI(3)**2
      DO 190 ILG=N+NP+MSTU(44)+MIN(ILC,10)+4,N+NP+MSTU(44)+5,-1
      IF(TDS.LE.P(ILG,4)) GOTO 200
      DO 190 J=1,4
  190 P(ILG+1,J)=P(ILG,J)
      ILG=N+NP+MSTU(44)+4
  200 DO 210 J=1,3
  210 P(ILG+1,J)=TDI(J)
      P(ILG+1,4)=TDS
  220 CONTINUE

C...Iterate direction of axis until stable maximum.
      P(N+NP+ILD,4)=0.
      ILG=0
  230 ILG=ILG+1
      THP=0.
  240 THPS=THP
      DO 250 J=1,3
      IF(THP.LE.1E-10) TDI(J)=P(N+NP+MSTU(44)+4+ILG,J)
      IF(THP.GT.1E-10) TDI(J)=TPR(J)
  250 TPR(J)=0.
      DO 260 I=N+1,N+NP
      SGN=SIGN(P(I,5),TDI(1)*P(I,1)+TDI(2)*P(I,2)+TDI(3)*P(I,3))
      DO 260 J=1,4-ILD
  260 TPR(J)=TPR(J)+SGN*P(I,J)
      THP=SQRT(TPR(1)**2+TPR(2)**2+TPR(3)**2)/PS
      IF(THP.GE.THPS+PARU(48)) GOTO 240

C...Save good axis. Try new initial axis until a number of tries agree.
      IF(THP.LT.P(N+NP+ILD,4)-PARU(48).AND.ILG.LT.MIN(10,NC)) GOTO 230
      IF(THP.GT.P(N+NP+ILD,4)+PARU(48)) THEN
        IAGR=0
        SGN=(-1.)**INT(RLU(0)+0.5)
        DO 270 J=1,3
  270   P(N+NP+ILD,J)=SGN*TPR(J)/(PS*THP)
        P(N+NP+ILD,4)=THP
        P(N+NP+ILD,5)=0.
      ENDIF
      IAGR=IAGR+1
  280 IF(IAGR.LT.MSTU(45).AND.ILG.LT.MIN(10,NC)) GOTO 230

C...Find minor axis and value by orthogonality.
      SGN=(-1.)**INT(RLU(0)+0.5)
      P(N+NP+3,1)=-SGN*P(N+NP+2,2)
      P(N+NP+3,2)=SGN*P(N+NP+2,1)
      P(N+NP+3,3)=0.
      THP=0.
      DO 290 I=N+1,N+NP
  290 THP=THP+P(I,5)*ABS(P(N+NP+3,1)*P(I,1)+P(N+NP+3,2)*P(I,2))
      P(N+NP+3,4)=THP/PS
      P(N+NP+3,5)=0.

C...Fill axis information. Rotate back to original coordinate system.
      DO 300 ILD=1,3
      K(N+ILD,1)=31
      K(N+ILD,2)=96
      K(N+ILD,3)=ILD
      K(N+ILD,4)=0
      K(N+ILD,5)=0
      DO 300 J=1,5
      P(N+ILD,J)=P(N+NP+ILD,J)
  300 V(N+ILD,J)=0.
      CALL LUDBRB(N+1,N+3,THE,PHI,0D0,0D0,0D0)

C...Calculate thrust and oblateness. Select storing option.
      THR=P(N+1,4)
      OBL=P(N+2,4)-P(N+3,4)
      MSTU(61)=N+1
      MSTU(62)=NP
      IF(MSTU(43).LE.1) MSTU(3)=3
      IF(MSTU(43).GE.2) N=N+3

      RETURN
      END
