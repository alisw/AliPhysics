 
C***********************************************************************
 
       SUBROUTINE PYTEST(MTEST)
 
C...Purpose: to provide a simple program (disguised as a subroutine) to
C...run at installation as a check that the program works as intended.
      COMMON/LUJETS/N,K(4000,5),P(4000,5),V(4000,5)
      COMMON/LUDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
      COMMON/LUDAT2/KCHG(500,3),PMAS(500,4),PARF(2000),VCKM(4,4)
      COMMON/LUDAT3/MDCY(500,3),MDME(2000,2),BRAT(2000),KFDP(2000,5)
      COMMON/PYSUBS/MSEL,MSUB(200),KFIN(2,-40:40),CKIN(200)
      COMMON/PYPARS/MSTP(200),PARP(200),MSTI(200),PARI(200)
      SAVE /LUJETS/,/LUDAT1/,/LUDAT2/,/LUDAT3/
      SAVE /PYSUBS/,/PYPARS/
 
C...Common initial values. Loop over initiating conditions.
      MSTP(122)=MAX(0,MIN(2,MTEST))
      MDCY(LUCOMP(111),1)=0
      NERR=0
      DO 130 IPROC=1,8
 
C...Reset process type, kinematics cuts, and the flags used.
      MSEL=0
      DO 100 ISUB=1,200
      MSUB(ISUB)=0
  100 CONTINUE
      CKIN(1)=2.
      CKIN(3)=0.
      MSTP(2)=1
      MSTP(11)=0
      MSTP(33)=0
      MSTP(81)=1
      MSTP(82)=1
      MSTP(111)=1
      MSTP(131)=0
      MSTP(133)=0
      PARP(131)=0.01
 
C...Prompt photon production at fixed target.
      IF(IPROC.EQ.1) THEN
        PZSUM=300.
        PESUM=SQRT(PZSUM**2+ULMASS(211)**2)+ULMASS(2212)
        PQSUM=2.
        MSEL=10
        CKIN(3)=5.
        CALL PYINIT('FIXT','pi+','p',PZSUM)
 
C...QCD processes at ISR energies.
      ELSEIF(IPROC.EQ.2) THEN
        PESUM=63.
        PZSUM=0.
        PQSUM=2.
        MSEL=1
        CKIN(3)=5.
        CALL PYINIT('CMS','p','p',PESUM)
 
C...W production + multiple interactions at CERN Collider.
      ELSEIF(IPROC.EQ.3) THEN
        PESUM=630.
        PZSUM=0.
        PQSUM=0.
        MSEL=12
        CKIN(1)=20.
        MSTP(82)=4
        MSTP(2)=2
        MSTP(33)=3
        CALL PYINIT('CMS','p','pbar',PESUM)
 
C...W/Z gauge boson pairs + pileup events at the Tevatron.
      ELSEIF(IPROC.EQ.4) THEN
        PESUM=1800.
        PZSUM=0.
        PQSUM=0.
        MSUB(22)=1
        MSUB(23)=1
        MSUB(25)=1
        CKIN(1)=200.
        MSTP(111)=0
        MSTP(131)=1
        MSTP(133)=2
        PARP(131)=0.04
        CALL PYINIT('CMS','p','pbar',PESUM)
 
C...Higgs production at LHC.
      ELSEIF(IPROC.EQ.5) THEN
        PESUM=15400.
        PZSUM=0.
        PQSUM=2.
        MSUB(3)=1
        MSUB(102)=1
        MSUB(123)=1
        MSUB(124)=1
        PMAS(25,1)=300.
        CKIN(1)=200.
        MSTP(81)=0
        MSTP(111)=0
        CALL PYINIT('CMS','p','p',PESUM)
 
C...Z' production at SSC.
      ELSEIF(IPROC.EQ.6) THEN
        PESUM=40000.
        PZSUM=0.
        PQSUM=2.
        MSEL=21
        PMAS(32,1)=600.
        CKIN(1)=400.
        MSTP(81)=0
        MSTP(111)=0
        CALL PYINIT('CMS','p','p',PESUM)
 
C...W pair production at 1 TeV e+e- collider.
      ELSEIF(IPROC.EQ.7) THEN
        PESUM=1000.
        PZSUM=0.
        PQSUM=0.
        MSUB(25)=1
        MSUB(69)=1
        MSTP(11)=1
        CALL PYINIT('CMS','e+','e-',PESUM)
 
C...Deep inelastic scattering at a LEP+LHC ep collider.
      ELSEIF(IPROC.EQ.8) THEN
        P(1,1)=0.
        P(1,2)=0.
        P(1,3)=8000.
        P(2,1)=0.
        P(2,2)=0.
        P(2,3)=-80.
        PESUM=8080.
        PZSUM=7920.
        PQSUM=0.
        MSUB(10)=1
        CKIN(3)=50.
        MSTP(111)=0
        CALL PYINIT('USER','p','e-',PESUM)
      ENDIF
 
C...Generate 20 events of each required type.
      DO 120 IEV=1,20
      CALL PYEVNT
      PESUMM=PESUM
      IF(IPROC.EQ.4) PESUMM=MSTI(41)*PESUM
 
C...Check conservation of energy/momentum/flavour.
      MERR=0
      DEVE=ABS(PLU(0,4)-PESUMM)+ABS(PLU(0,3)-PZSUM)
      DEVT=ABS(PLU(0,1))+ABS(PLU(0,2))
      DEVQ=ABS(PLU(0,6)-PQSUM)
      IF(DEVE.GT.2E-3*PESUM.OR.DEVT.GT.MAX(0.01,1E-4*PESUM).OR.
     &DEVQ.GT.0.1) MERR=1
      IF(MERR.NE.0) WRITE(MSTU(11),5000) IPROC,IEV
 
C...Check that all KF codes are known ones, and that partons/particles
C...satisfy energy-momentum-mass relation.
      DO 110 I=1,N
      IF(K(I,1).GT.20) GOTO 110
      IF(LUCOMP(K(I,2)).EQ.0) THEN
        WRITE(MSTU(11),5100) I
        MERR=MERR+1
      ENDIF
      PD=P(I,4)**2-P(I,1)**2-P(I,2)**2-P(I,3)**2-P(I,5)**2*
     &SIGN(1.,P(I,5))
      IF(ABS(PD).GT.MAX(0.1,0.002*P(I,4)**2,0.002*P(I,5)**2).OR.
     &(P(I,5).GE.0..AND.P(I,4).LT.0.)) THEN
        WRITE(MSTU(11),5200) I
        MERR=MERR+1
      ENDIF
  110 CONTINUE
 
C...Listing of erroneous events, and first event of each type.
      IF(MERR.GE.1) NERR=NERR+1
      IF(NERR.GE.10) THEN
        WRITE(MSTU(11),5300)
        CALL LULIST(1)
        STOP
      ENDIF
      IF(MTEST.GE.1.AND.(MERR.GE.1.OR.IEV.EQ.1)) THEN
        IF(MERR.GE.1) WRITE(MSTU(11),5400)
        CALL LULIST(1)
      ENDIF
  120 CONTINUE
 
C...List statistics for each process type.
      IF(MTEST.GE.1) CALL PYSTAT(1)
  130 CONTINUE
 
C...Summarize result of run.
      IF(NERR.EQ.0) WRITE(MSTU(11),5500)
      IF(NERR.GT.0) WRITE(MSTU(11),5600) NERR
      RETURN
 
C...Formats for information.
 5000 FORMAT(/5X,'Energy/momentum/flavour nonconservation for process',
     &I2,', event',I4)
 5100 FORMAT(/5X,'Entry no.',I4,' in following event not known code')
 5200 FORMAT(/5X,'Entry no.',I4,' in following event has faulty ',
     &'kinematics')
 5300 FORMAT(/5X,'This is the tenth error experienced! Something is ',
     &'wrong.'/5X,'Execution will be stopped after listing of event.')
 5400 FORMAT(5X,'Faulty event follows:')
 5500 FORMAT(//5X,'End result of PYTEST: no errors detected.')
 5600 FORMAT(//5X,'End result of PYTEST:',I2,' errors detected.'/
     &5X,'This should not have happened!')
      END
