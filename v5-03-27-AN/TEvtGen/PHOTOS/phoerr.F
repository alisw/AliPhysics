      SUBROUTINE PHOERR(IMES,TEXT,DATA)
C.----------------------------------------------------------------------
C.
C.    PHOTOS:   PHOton radiation in decays ERRror handling
C.
C.    Purpose:  Inform user  about (fatal) errors and warnings generated
C.              by either the user or the program.
C.
C.    Input Parameters:   IMES, TEXT, DATA
C.
C.    Output Parameters:  None
C.
C.    Author(s):  B. van Eijk                     Created at:  29/11/89
C.                                                Last Update: 10/01/92
C.
C.----------------------------------------------------------------------
      IMPLICIT NONE
      DOUBLE PRECISION DATA
      INTEGER IMES,IERROR
      REAL*8 SDATA
      INTEGER PHLUN
      COMMON/PHOLUN/PHLUN
      INTEGER PHOMES
      PARAMETER (PHOMES=10)
      INTEGER STATUS
      COMMON/PHOSTA/STATUS(PHOMES)
      CHARACTER TEXT*(*)
      SAVE IERROR
C--   security STOP switch  
      LOGICAL ISEC
      SAVE ISEC
      DATA ISEC /.TRUE./
      DATA IERROR/ 0/
      IF (IMES.LE.PHOMES) STATUS(IMES)=STATUS(IMES)+1
C--
C--   Count number of non-fatal errors...
      IF ((IMES.EQ. 6).AND.(STATUS(IMES).GE.2)) RETURN
      IF ((IMES.EQ.10).AND.(STATUS(IMES).GE.2)) RETURN
      SDATA=DATA
      WRITE(PHLUN,9000)
      WRITE(PHLUN,9120)
      GOTO (10,20,30,40,50,60,70,80,90,100),IMES
      WRITE(PHLUN,9130) IMES
      GOTO 120
   10 WRITE(PHLUN,9010) TEXT,INT(SDATA)
      GOTO 110
   20 WRITE(PHLUN,9020) TEXT,SDATA
      GOTO 110
   30 WRITE(PHLUN,9030) TEXT,SDATA
      GOTO 110
   40 WRITE(PHLUN,9040) TEXT
      GOTO 110
   50 WRITE(PHLUN,9050) TEXT,INT(SDATA)
      GOTO 110
   60 WRITE(PHLUN,9060) TEXT,SDATA
      GOTO 130
   70 WRITE(PHLUN,9070) TEXT,INT(SDATA)
      GOTO 110
   80 WRITE(PHLUN,9080) TEXT,INT(SDATA)
      GOTO 110
   90 WRITE(PHLUN,9090) TEXT,INT(SDATA)
      GOTO 110
  100 WRITE(PHLUN,9100) TEXT,SDATA
      GOTO 130
  110 CONTINUE
      WRITE(PHLUN,9140)
      WRITE(PHLUN,9120)
      WRITE(PHLUN,9000)
      IF (ISEC) THEN 
        STOP
      ELSE
        GOTO 130
      ENDIF
  120 IERROR=IERROR+1
      IF (IERROR.GE.10) THEN
        WRITE(PHLUN,9150)
        WRITE(PHLUN,9120)
        WRITE(PHLUN,9000)
        IF (ISEC) THEN 
          STOP
        ELSE
          GOTO 130
        ENDIF
      ENDIF
  130 WRITE(PHLUN,9120)
      WRITE(PHLUN,9000)
      RETURN
 9000 FORMAT(1H ,80('*'))
 9010 FORMAT(1H ,'* ',A,': Too many charged Particles, NCHARG =',I6,T81,
     &'*')
 9020 FORMAT(1H ,'* ',A,': Too much Bremsstrahlung required, PRSOFT = ',
     &F15.6,T81,'*')
 9030 FORMAT(1H ,'* ',A,': Combined Weight is exceeding 1., Weight = ',
     &F15.6,T81,'*')
 9040 FORMAT(1H ,'* ',A,
     &': Error in Rescaling charged and neutral Vectors',T81,'*')
 9050 FORMAT(1H ,'* ',A,
     &': Non matching charged Particle Pointer, NCHARG = ',I5,T81,'*')
 9060 FORMAT(1H ,'* ',A,
     &': Do you really work with a Particle of Spin: ',F4.1,' ?',T81,
     &'*')
 9070 FORMAT(1H ,'* ',A, ': Stack Length exceeded, NSTACK = ',I5 ,T81,
     &'*')
 9080 FORMAT(1H ,'* ',A,
     &': Random Number Generator Seed(1) out of Range: ',I8,T81,'*')
 9090 FORMAT(1H ,'* ',A,
     &': Random Number Generator Seed(2) out of Range: ',I8,T81,'*')
 9100 FORMAT(1H ,'* ',A,
     &': Available Phase Space below Cut-off: ',F15.6,' GeV/c^2',T81,
     &'*')
 9120 FORMAT(1H ,'*',T81,'*')
 9130 FORMAT(1H ,'* Funny Error Message: ',I4,' ! What to do ?',T81,'*')
 9140 FORMAT(1H ,'* Fatal Error Message, I stop this Run !',T81,'*')
 9150 FORMAT(1H ,'* 10 Error Messages generated, I stop this Run !',T81,
     &'*')
      END
