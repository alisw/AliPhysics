 
C*********************************************************************
 
      SUBROUTINE PYUPIN(ISUB,TITLE,SIGMAX)
 
C...Routine to be called by user to set up a user-defined process.
      COMMON/LUDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
      COMMON/PYINT2/ISET(200),KFPR(200,2),COEF(200,20),ICOL(40,4,2)
      COMMON/PYINT6/PROC(0:200)
      CHARACTER PROC*28
      SAVE /LUDAT1/,/PYINT2/,/PYINT6/
      CHARACTER*(*) TITLE
 
C...Check that subprocess number free.
      IF(ISUB.LT.1.OR.ISUB.GT.200.OR.ISET(ISUB).GE.0) THEN
        WRITE(MSTU(11),5000) ISUB
        STOP
      ENDIF
 
C...Fill information on new process.
      ISET(ISUB)=11
      COEF(ISUB,1)=SIGMAX
      PROC(ISUB)=TITLE//' '
 
C...Format for error output.
 5000 FORMAT(1X,'Error: user-defined subprocess code ',I4,
     &' not allowed.'//1X,'Execution stopped!')
 
      RETURN
      END
