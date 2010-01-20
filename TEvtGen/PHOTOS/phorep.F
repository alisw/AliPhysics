      SUBROUTINE PHOREP
C.----------------------------------------------------------------------
C.
C.    PHOTOS:   PHOton radiation in decays run summary REPort
C.
C.    Purpose:  Inform user about success and/or restrictions of PHOTOS
C.              encountered during execution.
C.
C.    Input Parameters:   Common /PHOSTA/
C.
C.    Output Parameters:  None
C.
C.    Author(s):  B. van Eijk                     Created at:  10/01/92
C.                                                Last Update: 10/01/92
C.
C.----------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER PHLUN
      COMMON/PHOLUN/PHLUN
      INTEGER PHOMES
      PARAMETER (PHOMES=10)
      INTEGER STATUS
      COMMON/PHOSTA/STATUS(PHOMES)
      INTEGER I
      LOGICAL ERROR
      ERROR=.FALSE.
      WRITE(PHLUN,9000)
      WRITE(PHLUN,9010)
      WRITE(PHLUN,9020)
      WRITE(PHLUN,9030)
      WRITE(PHLUN,9040)
      WRITE(PHLUN,9030)
      WRITE(PHLUN,9020)
      DO 10 I=1,PHOMES
        IF (STATUS(I).EQ.0) GOTO 10
        IF ((I.EQ.6).OR.(I.EQ.10)) THEN
          WRITE(PHLUN,9050) I,STATUS(I)
        ELSE
          ERROR=.TRUE.
          WRITE(PHLUN,9060) I,STATUS(I)
        ENDIF
   10 CONTINUE
      IF (.NOT.ERROR) WRITE(PHLUN,9070)
      WRITE(PHLUN,9020)
      WRITE(PHLUN,9010)
      RETURN
 9000 FORMAT(1H1)
 9010 FORMAT(1H ,80('*'))
 9020 FORMAT(1H ,'*',T81,'*')
 9030 FORMAT(1H ,'*',26X,25('='),T81,'*')
 9040 FORMAT(1H ,'*',30X,'PHOTOS Run Summary',T81,'*')
 9050 FORMAT(1H ,'*',22X,'Warning #',I2,' occured',I6,' times',T81,'*')
 9060 FORMAT(1H ,'*',23X,'Error #',I2,' occured',I6,' times',T81,'*')
 9070 FORMAT(1H ,'*',16X,'PHOTOS Execution has successfully terminated',
     &T81,'*')
      END
