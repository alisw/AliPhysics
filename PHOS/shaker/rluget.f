*CMZ :          17/07/98  15.44.33  by  Federico Carminati
*-- Author :
C*********************************************************************

      SUBROUTINE RLUGET(LFN,MOVE)

C...Purpose: to dump the state of the random number generator on a file
C...for subsequent startup from this state onwards.
*KEEP,LUDATR.
      COMMON /LUDATR/ MRLU(6),RRLU(100)
      SAVE /LUDATR/
*KEND.
      CHARACTER CHERR*8

C...Backspace required number of records (or as many as there are).
      IF(MOVE.LT.0) THEN
        NBCK=MIN(MRLU(6),-MOVE)
        DO 100 IBCK=1,NBCK
  100   BACKSPACE(LFN,ERR=110,IOSTAT=IERR)
        MRLU(6)=MRLU(6)-NBCK
      ENDIF

C...Unformatted write on unit LFN.
      WRITE(LFN,ERR=110,IOSTAT=IERR) (MRLU(I1),I1=1,5),
     &(RRLU(I2),I2=1,100)
      MRLU(6)=MRLU(6)+1
      RETURN

C...Write error.
  110 WRITE(CHERR,'(I8)') IERR
      CALL LUERRM(18,'(RLUGET:) error when accessing file, IOSTAT ='//
     &CHERR)

      RETURN
      END
