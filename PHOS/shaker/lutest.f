*CMZ :          17/07/98  15.44.36  by  Federico Carminati
*-- Author :
C*********************************************************************

      SUBROUTINE LUTEST(MTEST)

C...Purpose: to provide a simple program (disguised as subroutine) to
C...run at installation as a check that the program works as intended.
*KEEP,LUJETS.
      COMMON /LUJETS/ N,K(200000,5),P(200000,5),V(200000,5)
      SAVE /LUJETS/
*KEEP,LUDAT1.
      COMMON /LUDAT1/ MSTU(200),PARU(200),MSTJ(200),PARJ(200)
      SAVE /LUDAT1/
*KEND.
      DIMENSION PSUM(5),PINI(6),PFIN(6)

C...Loop over events to be generated.
      IF(MTEST.GE.1) CALL LUTABU(20)
      NERR=0
      DO 170 IEV=1,600

C...Reset parameter values. Switch on some nonstandard features.
      MSTJ(1)=1
      MSTJ(3)=0
      MSTJ(11)=1
      MSTJ(42)=2
      MSTJ(43)=4
      MSTJ(44)=2
      PARJ(17)=0.1
      PARJ(22)=1.5
      PARJ(43)=1.
      PARJ(54)=-0.05
      MSTJ(101)=5
      MSTJ(104)=5
      MSTJ(105)=0
      MSTJ(107)=1
      IF(IEV.EQ.301.OR.IEV.EQ.351.OR.IEV.EQ.401) MSTJ(116)=3

C...Ten events each for some single jets configurations.
      IF(IEV.LE.50) THEN
        ITY=(IEV+9)/10
        MSTJ(3)=-1
        IF(ITY.EQ.3.OR.ITY.EQ.4) MSTJ(11)=2
        IF(ITY.EQ.1) CALL LU1ENT(1,1,15.,0.,0.)
        IF(ITY.EQ.2) CALL LU1ENT(1,3101,15.,0.,0.)
        IF(ITY.EQ.3) CALL LU1ENT(1,-2203,15.,0.,0.)
        IF(ITY.EQ.4) CALL LU1ENT(1,-4,30.,0.,0.)
        IF(ITY.EQ.5) CALL LU1ENT(1,21,15.,0.,0.)

C...Ten events each for some simple jet systems; string fragmentation.
      ELSEIF(IEV.LE.130) THEN
        ITY=(IEV-41)/10
        IF(ITY.EQ.1) CALL LU2ENT(1,1,-1,40.)
        IF(ITY.EQ.2) CALL LU2ENT(1,4,-4,30.)
        IF(ITY.EQ.3) CALL LU2ENT(1,2,2103,100.)
        IF(ITY.EQ.4) CALL LU2ENT(1,21,21,40.)
        IF(ITY.EQ.5) CALL LU3ENT(1,2101,21,-3203,30.,0.6,0.8)
        IF(ITY.EQ.6) CALL LU3ENT(1,5,21,-5,40.,0.9,0.8)
        IF(ITY.EQ.7) CALL LU3ENT(1,21,21,21,60.,0.7,0.5)
        IF(ITY.EQ.8) CALL LU4ENT(1,2,21,21,-2,40.,0.4,0.64,0.6,0.12,0.2)

C...Seventy events with independent fragmentation and momentum cons.
      ELSEIF(IEV.LE.200) THEN
        ITY=1+(IEV-131)/16
        MSTJ(2)=1+MOD(IEV-131,4)
        MSTJ(3)=1+MOD((IEV-131)/4,4)
        IF(ITY.EQ.1) CALL LU2ENT(1,4,-5,40.)
        IF(ITY.EQ.2) CALL LU3ENT(1,3,21,-3,40.,0.9,0.4)
        IF(ITY.EQ.3) CALL LU4ENT(1,2,21,21,-2,40.,0.4,0.64,0.6,0.12,0.2)
        IF(ITY.GE.4) CALL LU4ENT(1,2,-3,3,-2,40.,0.4,0.64,0.6,0.12,0.2)

C...A hundred events with random jets (check invariant mass).
      ELSEIF(IEV.LE.300) THEN
  100   DO 110 J=1,5
  110   PSUM(J)=0.
        NJET=2.+6.*RLU(0)
        DO 120 I=1,NJET
        KFL=21
        IF(I.EQ.1) KFL=INT(1.+4.*RLU(0))
        IF(I.EQ.NJET) KFL=-INT(1.+4.*RLU(0))
        EJET=5.+20.*RLU(0)
        THETA=ACOS(2.*RLU(0)-1.)
        PHI=6.2832*RLU(0)
        IF(I.LT.NJET) CALL LU1ENT(-I,KFL,EJET,THETA,PHI)
        IF(I.EQ.NJET) CALL LU1ENT(I,KFL,EJET,THETA,PHI)
        IF(I.EQ.1.OR.I.EQ.NJET) PSUM(5)=PSUM(5)+ULMASS(KFL)
        DO 120 J=1,4
  120   PSUM(J)=PSUM(J)+P(I,J)
        IF(PSUM(4)**2-PSUM(1)**2-PSUM(2)**2-PSUM(3)**2.LT.
     &  (PSUM(5)+PARJ(32))**2) GOTO 100

C...Fifty e+e- continuum events with matrix elements.
      ELSEIF(IEV.LE.350) THEN
        MSTJ(101)=2
        CALL LUEEVT(0,40.)

C...Fifty e+e- continuum event with varying shower options.
      ELSEIF(IEV.LE.400) THEN
        MSTJ(42)=1+MOD(IEV,2)
        MSTJ(43)=1+MOD(IEV/2,4)
        MSTJ(44)=MOD(IEV/8,3)
        CALL LUEEVT(0,90.)

C...Fifty e+e- continuum events with coherent shower, including top.
      ELSEIF(IEV.LE.450) THEN
        MSTJ(104)=6
        CALL LUEEVT(0,500.)

C...Fifty Upsilon decays to ggg or gammagg with coherent shower.
      ELSEIF(IEV.LE.500) THEN
        CALL LUONIA(5,9.46)

C...One decay each for some heavy mesons.
      ELSEIF(IEV.LE.560) THEN
        ITY=IEV-501
        KFLS=2*(ITY/20)+1
        KFLB=8-MOD(ITY/5,4)
        KFLC=KFLB-MOD(ITY,5)
        CALL LU1ENT(1,100*KFLB+10*KFLC+KFLS,0.,0.,0.)

C...One decay each for some heavy baryons.
      ELSEIF(IEV.LE.600) THEN
        ITY=IEV-561
        KFLS=2*(ITY/20)+2
        KFLA=8-MOD(ITY/5,4)
        KFLB=KFLA-MOD(ITY,5)
        KFLC=MAX(1,KFLB-1)
        CALL LU1ENT(1,1000*KFLA+100*KFLB+10*KFLC+KFLS,0.,0.,0.)
      ENDIF

C...Generate event. Find total momentum, energy and charge.
      DO 130 J=1,4
  130 PINI(J)=PLU(0,J)
      PINI(6)=PLU(0,6)
      CALL LUEXEC
      DO 140 J=1,4
  140 PFIN(J)=PLU(0,J)
      PFIN(6)=PLU(0,6)

C...Check conservation of energy, momentum and charge;
C...usually exact, but only approximate for single jets.
      MERR=0
      IF(IEV.LE.50) THEN
        IF((PFIN(1)-PINI(1))**2+(PFIN(2)-PINI(2))**2.GE.4.) MERR=MERR+1
        EPZREM=PINI(4)+PINI(3)-PFIN(4)-PFIN(3)
        IF(EPZREM.LT.0..OR.EPZREM.GT.2.*PARJ(31)) MERR=MERR+1
        IF(ABS(PFIN(6)-PINI(6)).GT.2.1) MERR=MERR+1
      ELSE
        DO 150 J=1,4
  150   IF(ABS(PFIN(J)-PINI(J)).GT.0001*PINI(4)) MERR=MERR+1
        IF(ABS(PFIN(6)-PINI(6)).GT.0.1) MERR=MERR+1
      ENDIF
      IF(MERR.NE.0) WRITE(MSTU(11),1000) (PINI(J),J=1,4),PINI(6),
     &(PFIN(J),J=1,4),PFIN(6)

C...Check that all KF codes are known ones, and that partons/particles
C...satisfy energy-momentum-mass relation. Store particle statistics.
      DO 160 I=1,N
      IF(K(I,1).GT.20) GOTO 160
      IF(LUCOMP(K(I,2)).EQ.0) THEN
        WRITE(MSTU(11),1100) I
        MERR=MERR+1
      ENDIF
      PD=P(I,4)**2-P(I,1)**2-P(I,2)**2-P(I,3)**2-P(I,5)**2
      IF(ABS(PD).GT.MAX(0.1,0.001*P(I,4)**2).OR.P(I,4).LT.0.) THEN
        WRITE(MSTU(11),1200) I
        MERR=MERR+1
      ENDIF
  160 CONTINUE
      IF(MTEST.GE.1) CALL LUTABU(21)

C...List all erroneous events and some normal ones.
      IF(MERR.NE.0.OR.MSTU(24).NE.0.OR.MSTU(28).NE.0) THEN
        CALL LULIST(2)
      ELSEIF(MTEST.GE.1.AND.MOD(IEV-5,100).EQ.0) THEN
        CALL LULIST(1)
      ENDIF

C...Stop execution if too many errors. Endresult of run.
      IF(MERR.NE.0) NERR=NERR+1
      IF(NERR.GE.10) THEN
        WRITE(MSTU(11),1300) IEV
        STOP
      ENDIF
  170 CONTINUE
      IF(MTEST.GE.1) CALL LUTABU(22)
      WRITE(MSTU(11),1400) NERR

C...Reset commonblock variables changed during run.
      MSTJ(2)=3
      PARJ(17)=0.
      PARJ(22)=1.
      PARJ(43)=0.5
      PARJ(54)=0.
      MSTJ(105)=1
      MSTJ(107)=0

C...Format statements for output.
 1000 FORMAT(/' Momentum, energy and/or charge were not conserved ',
     &'in following event'/' sum of',9X,'px',11X,'py',11X,'pz',11X,
     &'E',8X,'charge'/' before',2X,4(1X,F12.5),1X,F8.2/' after',3X,
     &4(1X,F12.5),1X,F8.2)
 1100 FORMAT(/5X,'Entry no.',I4,' in following event not known code')
 1200 FORMAT(/5X,'Entry no.',I4,' in following event has faulty ',
     &'kinematics')
 1300 FORMAT(/5X,'Ten errors experienced by event ',I3/
     &5X,'Something is seriously wrong! Execution stopped now!')
 1400 FORMAT(/5X,'Number of erroneous or suspect events in run:',I3/
     &5X,'(0 fine, 1 acceptable if a single jet, ',
     &'>=2 something is wrong)')

      RETURN
      END
