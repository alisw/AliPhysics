 
C*********************************************************************
 
      SUBROUTINE PYUPEV(ISUB,SIGEV)
 
C...Dummy routine, to be replaced by user. When called from PYTHIA
C...the subprocess number ISUB will be given, and PYUPEV is supposed
C...to generate an event of this type, to be stored in the PYUPPR
C...commonblock. SIGEV gives the differential cross-section associated
C...with the event, i.e. the acceptance probability of the event is
C...taken to be SIGEV/SIGMAX, where SIGMAX was given in the PYUPIN
C...call.
      COMMON/LUDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
      SAVE /LUDAT1/
      COMMON/PYUPPR/NUP,KUP(20,7),PUP(20,5),NFUP,IFUP(10,2),Q2UP(0:10)
      SAVE /PYUPPR/
 
C...Stop program if this routine is ever called.
C...You should not copy these lines to your own routine.
      WRITE(MSTU(11),5000)
      IF(RLU(0).LT.10.) STOP
      SIGEV=ISUB
 
C...Format for error printout.
 5000 FORMAT(1X,'Error: you did not link your PYUPEV routine ',
     &'correctly.'/1X,'Dummy routine in PYTHIA file called instead.'/
     &1X,'Execution stopped!')
 
      RETURN
      END
