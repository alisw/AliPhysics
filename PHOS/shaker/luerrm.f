*CMZ :          17/07/98  15.44.33  by  Federico Carminati
*-- Author :
C*********************************************************************

      SUBROUTINE LUERRM(MERR,CHMESS)

C...Purpose: to inform user of errors in program execution.
*KEEP,LUJETS.
      COMMON /LUJETS/ N,K(200000,5),P(200000,5),V(200000,5)
      SAVE /LUJETS/
*KEEP,LUDAT1.
      COMMON /LUDAT1/ MSTU(200),PARU(200),MSTJ(200),PARJ(200)
      SAVE /LUDAT1/
*KEND.
      CHARACTER CHMESS*(*)

C...Write first few warnings, then be silent.
      IF(MERR.LE.10) THEN
        MSTU(27)=MSTU(27)+1
        MSTU(28)=MERR
        IF(MSTU(25).EQ.1.AND.MSTU(27).LE.MSTU(26)) WRITE(MSTU(11),1000)
     &  MERR,MSTU(31),CHMESS

C...Write first few errors, then be silent or stop program.
      ELSEIF(MERR.LE.20) THEN
        MSTU(23)=MSTU(23)+1
        MSTU(24)=MERR-10
        IF(MSTU(21).GE.1.AND.MSTU(23).LE.MSTU(22)) WRITE(MSTU(11),1100)
     &  MERR-10,MSTU(31),CHMESS
        IF(MSTU(21).GE.2.AND.MSTU(23).GT.MSTU(22)) THEN
          WRITE(MSTU(11),1100) MERR-10,MSTU(31),CHMESS
          WRITE(MSTU(11),1200)
          IF(MERR.NE.17) CALL LULIST(2)
          STOP
        ENDIF

C...Stop program in case of irreparable error.
      ELSE
        WRITE(MSTU(11),1300) MERR-20,MSTU(31),CHMESS
        STOP
      ENDIF

C...Formats for output.
 1000 FORMAT(/5X,'Advisory warning type',I2,' given after',I6,
     &' LUEXEC calls:'/5X,A)
 1100 FORMAT(/5X,'Error type',I2,' has occured after',I6,
     &' LUEXEC calls:'/5X,A)
 1200 FORMAT(5X,'Execution will be stopped after listing of last ',
     &'event!')
 1300 FORMAT(/5X,'Fatal error type',I2,' has occured after',I6,
     &' LUEXEC calls:'/5X,A/5X,'Execution will now be stopped!')

      RETURN
      END
