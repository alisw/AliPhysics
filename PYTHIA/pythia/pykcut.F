 
C*********************************************************************
 
      SUBROUTINE PYKCUT(MCUT)
 
C...Dummy routine, which the user can replace in order to make cuts on
C...the kinematics on the parton level before the matrix elements are
C...evaluated and the event is generated. The cross-section estimates
C...will automatically take these cuts into account, so the given
C...values are for the allowed phase space region only. MCUT=0 means
C...that the event has passed the cuts, MCUT=1 that it has failed.
      COMMON/LUDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
      SAVE /LUDAT1/
      COMMON/PYINT1/MINT(400),VINT(400)
      COMMON/PYINT2/ISET(200),KFPR(200,2),COEF(200,20),ICOL(40,4,2)
      SAVE /PYINT1/,/PYINT2/
 
C...Set default value (accepting event) for MCUT.
      MCUT=0
 
C...Read out subprocess number.
      ISUB=MINT(1)
      ISTSB=ISET(ISUB)
 
C...Read out tau, y*, cos(theta), tau' (where defined, else =0).
      TAU=VINT(21)
      YST=VINT(22)
      CTH=0.
      IF(ISTSB.EQ.2.OR.ISTSB.EQ.4.OR.ISTSB.EQ.6) CTH=VINT(23)
      TAUP=0.
      IF(ISTSB.GE.3.AND.ISTSB.LE.5) TAUP=VINT(26)
 
C...Calculate x_1, x_2, x_F.
      IF(ISTSB.LE.2.OR.ISTSB.GE.6) THEN
        X1=SQRT(TAU)*EXP(YST)
        X2=SQRT(TAU)*EXP(-YST)
      ELSE
        X1=SQRT(TAUP)*EXP(YST)
        X2=SQRT(TAUP)*EXP(-YST)
      ENDIF
      XF=X1-X2
 
C...Calculate shat, that, uhat, p_T^2.
      SHAT=TAU*VINT(2)
      SQM3=VINT(63)
      SQM4=VINT(64)
      RM3=SQM3/SHAT
      RM4=SQM4/SHAT
      BE34=SQRT(MAX(0.,(1.-RM3-RM4)**2-4.*RM3*RM4))
      RPTS=4.*VINT(71)**2/SHAT
      BE34L=SQRT(MAX(0.,(1.-RM3-RM4)**2-4.*RM3*RM4-RPTS))
      RM34=2.*RM3*RM4
      RSQM=1.+RM34
      RTHM=(4.*RM3*RM4+RPTS)/(1.-RM3-RM4+BE34L)
      THAT=-0.5*SHAT*MAX(RTHM,1.-RM3-RM4-BE34*CTH)
      UHAT=-0.5*SHAT*MAX(RTHM,1.-RM3-RM4+BE34*CTH)
      PT2=MAX(VINT(71)**2,0.25*SHAT*BE34**2*(1.-CTH**2))
 
C...Decisions by user to be put here.
 
C...Stop program if this routine is ever called.
C...You should not copy these lines to your own routine.
      WRITE(MSTU(11),5000)
      IF(RLU(0).LT.10.) STOP
 
C...Format for error printout.
 5000 FORMAT(1X,'Error: you did not link your PYKCUT routine ',
     &'correctly.'/1X,'Dummy routine in PYTHIA file called instead.'/
     &1X,'Execution stopped!')
 
      RETURN
      END
