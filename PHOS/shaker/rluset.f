*CMZ :          17/07/98  15.44.33  by  Federico Carminati
*-- Author :
C*********************************************************************

      SUBROUTINE RLUSET(LFN,MOVE)

C...Purpose: to read a state of the random number generator from a file
C...for subsequent generation from this state onwards.
*KEEP,LUDATR.
      COMMON /LUDATR/ MRLU(6),RRLU(100)
      SAVE /LUDATR/
*KEND.
      CHARACTER CHERR*8

C...Backspace required number of records (or as many as there are).
      IF(MOVE.LT.0) THEN
        NBCK=MIN(MRLU(6),-MOVE)
        DO 100 IBCK=1,NBCK
  100   BACKSPACE(LFN,ERR=120,IOSTAT=IERR)
        MRLU(6)=MRLU(6)-NBCK
      ENDIF

C...Unformatted read from unit LFN.
      NFOR=1+MAX(0,MOVE)
      DO 110 IFOR=1,NFOR
  110 READ(LFN,ERR=120,IOSTAT=IERR) (MRLU(I1),I1=1,5),
     &(RRLU(I2),I2=1,100)
      MRLU(6)=MRLU(6)+NFOR
      RETURN

C...Write error.
  120 WRITE(CHERR,'(I8)') IERR
      CALL LUERRM(18,'(RLUSET:) error when accessing file, IOSTAT ='//
     &CHERR)

      RETURN
      END
