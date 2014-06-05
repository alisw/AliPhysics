      DOUBLE PRECISION FUNCTION GFACT(QQ)
      IMPLICIT NONE
C factor G to be used as inteligent retabulation as in paper 
C Kuhn Santamaria
      COMMON / PARMAS / AMTAU,AMNUTA,AMEL,AMNUE,AMMU,AMNUMU
     *                 ,AMPIZ,AMPI,AMRO,GAMRO,AMA1,GAMA1
     *                 ,AMK,AMKZ,AMKST,GAMKST
C
      REAL*4            AMTAU,AMNUTA,AMEL,AMNUE,AMMU,AMNUMU
     *                 ,AMPIZ,AMPI,AMRO,GAMRO,AMA1,GAMA1
     *                 ,AMK,AMKZ,AMKST,GAMKST
      DOUBLE PRECISION DGAMQQ,QQ,THR,THU,TH0,LAM
      INTEGER II,INIT,INUM
      DOUBLE PRECISION AF4,AF5,AA1,AA2,AA4
      DOUBLE PRECISION F1,F2,F3,F4,F5
      DOUBLE PRECISION AA,BB,CC,DD,A,B,C,D,X
      SAVE A,B,C,D,AF4,AF5,AA,BB,CC,DD,THR,THU,TH0
      
      CALL IFGFACT(1,II,INIT)
      IF(INIT.EQ.0) THEN
       THR=(AMPI+AMRO)**2
       TH0=9*AMPI**2
       THU=THR
       CALL getFF3PISCAL(INUM)  ! we switch off part of the contr to 3 pi mode
       CALL setFF3PISCAL(0)
C below THR , we calculate at distance 1/4 1/2 and 1 between 
C minimum and maximum for this range
       LAM=(THR-TH0)/4
       AA1=DGAMQQ(TH0+  LAM)
       AA2=DGAMQQ(TH0+2*LAM)
       AA4=DGAMQQ(TH0+4*LAM)
C above THR we calculate at THR times 1 2 3 and 4 for higher range
       F1=DGAMQQ(THU)
       F2=DGAMQQ(1.5*THU)
       F3=DGAMQQ(2.*THU)
       F4=DGAMQQ(3*THU)
       F5=DGAMQQ(3.5*THU)
       CALL setFF3PISCAL(INUM) ! we switch back part of the contr to 3 pi mode

C we calculate coefs for expansion  ( polynomial of order 2 once 
C                                     X**2 factorized out, X=QQ-TH0 )

       AA1 = AA1/LAM**2            ! RESCALING due to factorized X**2
       AA2 = AA2/(2.*LAM)**2       !
       AA4 = AA4/(4.*LAM)**2       !

       BB=1./8.*(10.*AA2-AA4-16.*AA1)
       AA=(8.*AA1-AA2-4.*BB)/6.
       CC=AA1-AA-BB
       AA=AA/LAM
       BB=BB/LAM**2
       CC=CC/LAM**3
C calulate coefs, assuming it is polynomial order 3 note negative pwrs
       D=-9*(-4*F3-F1+4*F2+F4)
       C=3*(F4+F1-2*F3-11./18.*D)
       A=F3-F1+0.5*C+0.75*D
       B=F1-A-C-D
       A=A/THU
       B=B
       C=C*THU
       D=D*THU**2
       AF4=F4
       AF5=F5
c       write(*,*) "A=",AA,"B=",-BB,"C=",CC,"D=",A,"E=",-B,"F=",C
c       write(*,*) "G=",-D,"H=",AF4,"P=",AF5-AF4
      ENDIF
      IF (QQ.GT.3*THU) THEN
        GFACT=AF4+(AF5-AF4)*(QQ-3*THU)*2/THU
      ELSEIF (QQ.GT.THR) THEN
        GFACT=A*QQ+B+C/QQ+D/QQ**2
      ELSEIF(QQ.LE.TH0) THEN
        GFACT=0.0
      ELSE
        X=QQ-TH0
        GFACT=AA*X+BB*X**2+CC*X**3
        GFACT=X**2*GFACT
      ENDIF
      END


      DOUBLE PRECISION FUNCTION DGAMQQ(XQQB)
C **************************************************************
C     calculates \tau^- -> pi^- pi^- pi^+ nu width as function of QQ (XQQB)
C     formulas (19) of ref [2a] integration over S1
C     limit of integration (21) of ref [2a]  see also [4]
C     called from main function
C **************************************************************
      IMPLICIT NONE 
      COMMON/PRECINT/          EPSSQ,ABS1
      DOUBLE PRECISION         EPSSQ,ABS1
      COMMON /EXTERNAL/ XQQA
      DOUBLE PRECISION  XQQA
      EXTERNAL          GAUS,DGAMQQS1
      DOUBLE PRECISION  GAUS,DGAMQQS1
      DOUBLE PRECISION  XQQB,EPS,UPS1,DOWNS1
      COMMON / PARMAS / AMTAU,AMNUTA,AMEL,AMNUE,AMMU,AMNUMU
     *                 ,AMPIZ,AMPI,AMRO,GAMRO,AMA1,GAMA1
     *                 ,AMK,AMKZ,AMKST,GAMKST
C
      REAL*4            AMTAU,AMNUTA,AMEL,AMNUE,AMMU,AMNUMU
     *                 ,AMPIZ,AMPI,AMRO,GAMRO,AMA1,GAMA1
     *                 ,AMK,AMKZ,AMKST,GAMKST

      XQQA = XQQB
      EPS = EPSSQ/3.D0
      UPS1=(DSQRT(XQQB) - AMPI)**2-ABS1              ! limits on S1
      DOWNS1=4.D0*AMPI**2+ABS1

      DGAMQQ = GAUS(DGAMQQS1,DOWNS1,UPS1,EPS)

      RETURN
      END



      DOUBLE PRECISION FUNCTION DGAMQQS1(S1)
      IMPLICIT NONE
      DOUBLE PRECISION                   S1
C **************************************************************
C     calculates \tau^- -> pi^- pi^- pi^+ nu width 
C           as function of QQ,S1 
C     GAUS integrant in DGAMQQ (XQQA- hidden argument)
C     calculates tau^- -> pi^- pi^- pi^+ nu spectrum as function of S1
C     formulas (19) of ref [2a] see also [4]
C     limit of integration (21) of ref [2a]
C **************************************************************
      EXTERNAL          GAUS2,FFWID3PI,DGAMQQS1S3
      DOUBLE PRECISION  GAUS2,FFWID3PI,DGAMQQS1S3
      COMMON /INTERNAL/ S1A
      DOUBLE PRECISION  S1A
      COMMON /EXTERNAL/   XQQA
      DOUBLE PRECISION    XQQA
      COMMON/PRECINT/          EPSSQ,ABS1
      DOUBLE PRECISION         EPSSQ,ABS1
      DOUBLE PRECISION EPS,UPS3,DOWNS3
      DOUBLE PRECISION XLAM,X,Y,Z
      DOUBLE PRECISION XAMPI2
      COMMON / PARMAS / AMTAU,AMNUTA,AMEL,AMNUE,AMMU,AMNUMU
     *                 ,AMPIZ,AMPI,AMRO,GAMRO,AMA1,GAMA1
     *                 ,AMK,AMKZ,AMKST,GAMKST
C
      REAL*4            AMTAU,AMNUTA,AMEL,AMNUE,AMMU,AMNUMU
     *                 ,AMPIZ,AMPI,AMRO,GAMRO,AMA1,GAMA1
     *                 ,AMK,AMKZ,AMKST,GAMKST

      Xlam(x,y,z) = sqrt(abs((x-y-z)**2 - 4.*y*z))
      S1A = S1
      EPS = EPSSQ/9.D0
      XAMPI2=AMPI**2
      
      UPS3 = (XQQA - AMPI**2)**2 -                              ! limits on S3
     &    ( XLAM(XQQA,S1,XAMPI2) 
     &         - XLAM(S1,XAMPI2,XAMPI2) )**2
      DOWNS3 = (XQQA - AMPI**2)**2 - 
     &    (XLAM(XQQA,S1,XAMPI2) 
     &          + XLAM(S1,XAMPI2,XAMPI2) )**2

      UPS3 = UPS3/4./S1
      DOWNS3 = DOWNS3/4./S1
      
      DGAMQQS1 = GAUS2(DGAMQQS1S3,DOWNS3,UPS3,EPS)

      RETURN
      END


      DOUBLE PRECISION FUNCTION DGAMQQS1S3(XS3)
      IMPLICIT NONE
C **************************************************************
C     calculates \tau^- -> pi^- pi^- pi^+ nu width 
C           as function of QQ,S1,S3
C **************************************************************
C      EXTERNAL          GAUS2,FFWID3PI
      DOUBLE PRECISION  FFWID3PI,XS3
      COMMON /INTERNAL/ XS1A
      DOUBLE PRECISION  XS1A
      COMMON /EXTERNAL/   XQQA
      DOUBLE PRECISION    XQQA
      COMMON / PARMAS / AMTAU,AMNUTA,AMEL,AMNUE,AMMU,AMNUMU
     *                 ,AMPIZ,AMPI,AMRO,GAMRO,AMA1,GAMA1
     *                 ,AMK,AMKZ,AMKST,GAMKST
C
      REAL*4            AMTAU,AMNUTA,AMEL,AMNUE,AMMU,AMNUMU
     *                 ,AMPIZ,AMPI,AMRO,GAMRO,AMA1,GAMA1
     *                 ,AMK,AMKZ,AMKST,GAMKST


      DGAMQQS1S3 = FFWID3PI(XQQA,XS1A,XS3)

      RETURN

      END
