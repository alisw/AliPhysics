 
C*********************************************************************
 
      SUBROUTINE PDFSET(PARM,VALUE)
 
C...Dummy routine, to be removed when PDFLIB is to be linked.
      COMMON/LUDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
      SAVE /LUDAT1/
      CHARACTER*20 PARM(20)
      DOUBLE PRECISION VALUE(20)
 
C...Stop program if this routine is ever called.
      WRITE(MSTU(11),5000)
      IF(RLU(0).LT.10.) STOP
      PARM(20)=PARM(1)
      VALUE(20)=VALUE(1)
 
C...Format for error printout.
 5000 FORMAT(1X,'Error: you did not link PDFLIB correctly.'/
     &1X,'Dummy routine PDFSET in PYTHIA file called instead.'/
     &1X,'Execution stopped!')
 
      RETURN
      END
