      DOUBLE PRECISION FUNCTION FFWID3PI(QQ,S1,S3)
      IMPLICIT NONE      
      DOUBLE PRECISION                   QQ,S1,S3
C **************************************************************
C     Input:  QQ S1 S3   ! mpi-pi-pi+**2   mpi-pi+**2  mpi-pi-**2
C     Calls: functions FORM1, FORM2, FORM4
C     Uses constants: tau mass, pi mass, normalization constant.
C     Load: Initialized tauola library, 
C     Output:  d\Gamma(tau --> 3pi nu)/(dQQ dS1 dS3)
C     Remark: If QQ S1 S3 are outside of the phase space 
C             function FFWID3PI returns zero.
C **************************************************************
      COMPLEX F1,F2,F4, FORM1,FORM2, FORM4
      DOUBLE PRECISION V11,V12,V22,GGF2,VUD2,ABS1,QQMIN,
     &                 QQMAX,S3MAX,S3MIN,S1MIN,S1MAX
      REAL XQQ,XS1,XS3, XS2,RQQ 
      COMMON / DECPAR / GFERMI,GV,GA,CCABIB,SCABIB,GAMEL
      REAL*4            GFERMI,GV,GA,CCABIB,SCABIB,GAMEL
      COMMON / PARMAS / AMTAU,AMNUTA,AMEL,AMNUE,AMMU,AMNUMU
     *                 ,AMPIZ,AMPI,AMRO,GAMRO,AMA1,GAMA1
     *                 ,AMK,AMKZ,AMKST,GAMKST
      DOUBLE PRECISION  XLAM,X,Y,Z
      DOUBLE PRECISION  XAMPI2
      DOUBLE PRECISION  GETFPIRPT
      REAL*4            AMTAU,AMNUTA,AMEL,AMNUE,AMMU,AMNUMU
     *                 ,AMPIZ,AMPI,AMRO,GAMRO,AMA1,GAMA1
     *                 ,AMK,AMKZ,AMKST,GAMKST
      DOUBLE PRECISION        PI   
      DATA                    PI /3.141592653589793238462643D0/
      INTEGER  IMODE,IDUM,IFRCHL
      REAL RRQ,RCHLWIDA1PI


      XLAM(X,Y,Z)= sqrt(abs((x-y-z)**2 - 4.*y*z))
     
      ABS1 = 1.d-5


      GGF2 = GFERMI**2
      VUD2 = CCABIB**2


C     TO CHANGE THE VARIABLES TO SINGLE PRECISION
C     INPUT FOR FORM1,FORM2,FORM4 IS SINGLE PRECISION

      XQQ = QQ
      XS1 = S1
      XS2 = QQ -S1-S3 + 3.*AMPI**2
      XS3 = S3
      XAMPI2 = AMPI**2
C     Limits for PHASE SPACE
C     Limits for QQ
      QQMIN = 9.D0*AMPI**2
      QQMAX  = (AMTAU-AMNUTA)**2

C     Limits for S1
      S1MAX=(DSQRT(QQ) - AMPI)**2 -ABS1
      S1MIN=4.D0*AMPI**2 +ABS1

C    LIMIT FOR XS3
      S3MAX = (QQ - AMPI**2)**2 - 
     &    ( XLAM(QQ,S1,XAMPI2) 
     &         - XLAM(S1,XAMPI2,XAMPI2) )**2
      S3MIN = (QQ - AMPI**2)**2 - 
     &    (XLAM(QQ,S1,XAMPI2) 
     &          + XLAM(S1,XAMPI2,XAMPI2) )**2

      S3MAX = S3MAX/4./S1
      S3MIN = S3MIN/4./S1

C     Check on PHASE SPACE
C 
      IF((XS2.LE.0.) .OR.(S3MAX.LE.S3MIN)
     &   .OR.(XS1.LE.S1MIN).OR.(XS1.GE.S1MAX)
     &   .OR.(XS3.LE.S3MIN).OR.(XS3.GE.S3MAX)
     &   .OR.(QQ.LE.QQMIN).OR.(QQ.GE.QQMAX)
     &    )  THEN 
        FFWID3PI = 0.D0
        RETURN
      ENDIF

C
      V11 = -XS1+4.D0*AMPI**2 -(XS2-XS3)**2/(4.D0*XQQ)
      V22 = -XS2+4.D0*AMPI**2 - (XS3-XS1)**2/(4.D0*XQQ)
      V12 = 0.5D0*(XS3-XS1-XS2+4.D0*AMPI**2)-0.25D0*(XS3-XS2)*(XS3-XS1)/XQQ


      F1 = FORM1(0,XQQ,XS1,XS2)
      F2 = FORM2(0,XQQ,XS2,XS1)
      F4 = FORM4(0,XQQ,XS2,XS1,XS3)

C formula 3.46 of [3]   
      FFWID3PI = ABS(F1*CONJG(F1))*V11+ABS(F2*CONJG(F2))*V22+
     $       2.D0*REAL(F1*CONJG(F2))*V12 

      CALL IFGFACT(2,IMODE,IDUM)
      CALL INIRChLget(IFRCHL)
      IF (IMODE.EQ.0) THEN
C      VERSION A: The 3 pion contribution to the a1 width  
C      factor of a1 phase space and zeroing a1 propagator etc.is 
C      done in RCHLWIDA1PI 
       RQQ=QQ
       IF (IFRCHL.EQ.1) THEN
C       CASE OF RCHL 
        FFWID3PI   = RCHLWIDA1PI(RQQ,FFWID3PI)
       ELSE
C       NOT READY YET
        WRITE(*,*) 'FFWID3PI is not ready for non rchl currents'
        STOP
       ENDIF
C      to get the total a1 width the contribution from (KKPI)- and K-K0pi0
C      channels have to be added
      ELSE
C      VERSION B: calculation of 3 pion spectra in tau to 3pi nu channel.

C      factor for phase space of tau to XQQ nu decay and contribution from F4
C      (formula 3.21 of [3]) 

      FFWID3PI = (- FFWID3PI/3.D0*(1.D0+2.D0*XQQ/(AMTAU**2))+ 
     $       XQQ*ABS(F4*CONJG(F4)))*(AMTAU**2/XQQ-1.D0)**2

C Flux factor and normalization const.
        FFWID3PI =     
     &       GGF2*VUD2/(128.D0*(2.D0*PI)**5*AMTAU)/2.d0
     &            *FFWID3PI  
       IF (IFRCHL.EQ.1) THEN
C       CASE OF RCHL 
C RChL normalization constant 
        FFWID3PI =FFWID3PI/GETFPIRPT(1)**2
       ELSE
C       NOT READY YET
        WRITE(*,*) 'FFWID3PI is not ready for non rchl currents'
        STOP
       ENDIF

      ENDIF

      RETURN
      END

      REAL FUNCTION RCHLWIDA1PI(RQQ,FFWID3PI)
C The 3 pion contribution to the a1 width   
C (in [3] simple pretabulation is used through formula 3.48)
C for calculation of g(QQ) of 3.45 3.46 of [3] in RChL style
C a1 propagator has to be taken with the zero width.

      IMPLICIT NONE
      COMPLEX FA1RCHL
      REAL RQQ
      DOUBLE PRECISION FFWID3PI
      COMMON / DECPAR / GFERMI,GV,GA,CCABIB,SCABIB,GAMEL
      REAL*4            GFERMI,GV,GA,CCABIB,SCABIB,GAMEL
      COMMON / PARMAS / AMTAU,AMNUTA,AMEL,AMNUE,AMMU,AMNUMU
     *                 ,AMPIZ,AMPI,AMRO,GAMRO,AMA1,GAMA1
     *                 ,AMK,AMKZ,AMKST,GAMKST
      DOUBLE PRECISION  XLAM,X,Y,Z
      DOUBLE PRECISION  XAMPI2
      DOUBLE PRECISION  GETFPIRPT
      REAL*4            AMTAU,AMNUTA,AMEL,AMNUE,AMMU,AMNUMU
     *                 ,AMPIZ,AMPI,AMRO,GAMRO,AMA1,GAMA1
     *                 ,AMK,AMKZ,AMKST,GAMKST

      COMMON/RCHT_3PI/        FPI_RPT,FV_RPT,GV_RPT,FA_RPT,BETA_RHO,FK_RPT
     &                       ,FV1_RPT,GV1_RPT
      DOUBLE PRECISION        FPI_RPT,FV_RPT,GV_RPT,FA_RPT,BETA_RHO,FK_RPT
     &                       ,FV1_RPT,GV1_RPT
      DOUBLE PRECISION        PI   
      DATA                    PI /3.141592653589793238462643D0/

C     AMA1 should be replaced by variable from the rchl namespace.
      RCHLWIDA1PI=- 1.0/REAL(FA1RCHL(RQQ)*CONJG(FA1RCHL(RQQ)))/RQQ**2 
     $            /(96.D0*8.D0*PI**3*AMA1)/(FA_RPT**2*FPI_RPT**2)
     $            *FFWID3PI/2.d0 
      END

      DOUBLE PRECISION FUNCTION  GETFPIRPT(I)
      IMPLICIT NONE
      COMMON/RCHT_3PI/        FPI_RPT,FV_RPT,GV_RPT,FA_RPT,BETA_RHO,FK_RPT
     &                       ,FV1_RPT,GV1_RPT
      DOUBLE PRECISION        FPI_RPT,FV_RPT,GV_RPT,FA_RPT,BETA_RHO,FK_RPT
     &                       ,FV1_RPT,GV1_RPT
      INTEGER I
      GETFPIRPT=FPI_RPT
      END
