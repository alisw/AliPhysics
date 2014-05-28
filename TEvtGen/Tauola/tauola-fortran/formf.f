

      FUNCTION FORMOM(XMAA,XMOM)
C     ==================================================================
C     formfactorfor pi-pi0 gamma final state
C      R. Decker, Z. Phys C36 (1987) 487.
C     ==================================================================
      COMMON / PARMAS / AMTAU,AMNUTA,AMEL,AMNUE,AMMU,AMNUMU
     *                 ,AMPIZ,AMPI,AMRO,GAMRO,AMA1,GAMA1
     *                 ,AMK,AMKZ,AMKST,GAMKST
C
      REAL*4            AMTAU,AMNUTA,AMEL,AMNUE,AMMU,AMNUMU
     *                 ,AMPIZ,AMPI,AMRO,GAMRO,AMA1,GAMA1
     *                 ,AMK,AMKZ,AMKST,GAMKST
      COMMON / DECPAR / GFERMI,GV,GA,CCABIB,SCABIB,GAMEL
      REAL*4            GFERMI,GV,GA,CCABIB,SCABIB,GAMEL
      COMMON /TESTA1/ KEYA1
      COMPLEX BWIGN,FORMOM
      DATA ICONT /1/
* THIS INLINE FUNCT. CALCULATES THE SCALAR PART OF THE PROPAGATOR
      BWIGN(XM,AM,GAMMA)=1./CMPLX(XM**2-AM**2,GAMMA*AM)
* HADRON CURRENT
      FRO  =0.266*AMRO**2
      ELPHA=- 0.1
      AMROP = 1.7
      GAMROP= 0.26
      AMOM  =0.782
      GAMOM =0.0085
      AROMEG= 1.0
      GCOUP=12.924
      GCOUP=GCOUP*AROMEG
      FQED  =SQRT(4.0*3.1415926535/137.03604)
      FORMOM=FQED*FRO**2/SQRT(2.0)*GCOUP**2*BWIGN(XMOM,AMOM,GAMOM)
     $     *(BWIGN(XMAA,AMRO,GAMRO)+ELPHA*BWIGN(XMAA,AMROP,GAMROP))
     $     *(BWIGN( 0.0,AMRO,GAMRO)+ELPHA*BWIGN( 0.0,AMROP,GAMROP))
      END



C=======================================================================
      COMPLEX FUNCTION FK1AB(XMSQ,INDX)
C     ==================================================================
C     complex form-factor for a1+a1prime.                       AJW 1/98
C     ==================================================================

      COMPLEX F1,F2,AMPA,AMPB
      INTEGER IFIRST,INDX
      DATA IFIRST/0/

      IF (IFIRST.EQ.0) THEN
        IFIRST = 1
        XM1 = PKORB(1,19)
        XG1 = PKORB(2,19)
        XM2 = PKORB(1,20)
        XG2 = PKORB(2,20)

        XM1SQ = XM1*XM1
        GF1 = GFUN(XM1SQ)
        GG1 = XM1*XG1/GF1
        XM2SQ = XM2*XM2
        GF2 = GFUN(XM2SQ)
        GG2 = XM2*XG2/GF2
      END IF

      IF (INDX.EQ.1) THEN
        AMPA = CMPLX(PKORB(3,81),0.)
        AMPB = CMPLX(PKORB(3,82),0.)
      ELSE IF (INDX.EQ.2) THEN
        AMPA = CMPLX(PKORB(3,83),0.)
        AMPB = CMPLX(PKORB(3,84),0.)
      ELSEIF (INDX.EQ.3) THEN
        AMPA = CMPLX(PKORB(3,85),0.)
        AMPB = CMPLX(PKORB(3,86),0.)
      ELSEIF (INDX.EQ.4) THEN
        AMPA = CMPLX(PKORB(3,87),0.)
        AMPB = CMPLX(PKORB(3,88),0.)
      END IF

      GF = GFUN(XMSQ)
      FG1 = GG1*GF
      FG2 = GG2*GF
      F1 = CMPLX(-XM1SQ,0.0)/CMPLX(XMSQ-XM1SQ,FG1)
      F2 = CMPLX(-XM2SQ,0.0)/CMPLX(XMSQ-XM2SQ,FG2)
      FK1AB = AMPA*F1+AMPB*F2

      RETURN
      END

      FUNCTION FORM1(MNUM,QQ,S1,SDWA)
C     ==================================================================
C     formfactorfor F1 for 3 scalar final state
C     R. Fisher, J. Wess and F. Wagner Z. Phys C3 (1980) 313
C     H. Georgi, Weak interactions and modern particle theory,
C     The Benjamin/Cummings Pub. Co., Inc. 1984.
C     R. Decker, E. Mirkes, R. Sauer, Z. Was Karlsruhe preprint TTP92-25
C     and erratum !!!!!!
C     ==================================================================
C
      COMPLEX FORM1,WIGNER,WIGFOR,FPIKM,BWIGM
      COMMON / PARMAS / AMTAU,AMNUTA,AMEL,AMNUE,AMMU,AMNUMU
     *                 ,AMPIZ,AMPI,AMRO,GAMRO,AMA1,GAMA1
     *                 ,AMK,AMKZ,AMKST,GAMKST
C
      REAL*4            AMTAU,AMNUTA,AMEL,AMNUE,AMMU,AMNUMU
     *                 ,AMPIZ,AMPI,AMRO,GAMRO,AMA1,GAMA1
     *                 ,AMK,AMKZ,AMKST,GAMKST
      COMMON /IPChT/ IVER
      INTEGER        IVER

      COMPLEX FORMA1,FORMK1,FORMRO,FORMKS
      COMPLEX FA1A1P,FK1AB,F3PI,F3PI_RCHT
C
      IF     (MNUM.EQ.0) THEN
C ------------  3 pi hadronic state (a1)
C       FORMRO = FPIKM(SQRT(S1),AMPI,AMPI)
C       FORMRO = F3PI(1,QQ,S1,SDWA)
C       FORMA1 = FA1A1P(QQ)
C       FORM1 = FORMA1*FORMRO
      IF (IVER.EQ.0) THEN
        FORM1 = F3PI(1,QQ,S1,SDWA)
      ELSE
        FORM1 = F3PI_RCHT(1,QQ,S1,SDWA)
      ENDIF

      ELSEIF (MNUM.EQ.1) THEN
C ------------ K- pi- K+ (K*0 K-)
       FORMKS = BWIGM(S1,AMKST,GAMKST,AMPI,AMK)
       FORMA1 = FA1A1P(QQ)
       FORM1 = FORMA1*FORMKS

      ELSEIF (MNUM.EQ.2) THEN
C ------------ K0 pi- K0B (K*- K0)
       FORMKS = BWIGM(S1,AMKST,GAMKST,AMPI,AMK)
       FORMA1 = FA1A1P(QQ)
       FORM1 = FORMA1*FORMKS

      ELSEIF (MNUM.EQ.3) THEN
C ------------ K- pi0 K0 (K*0 K-)
       FORMKS = BWIGM(S1,AMKST,GAMKST,AMPI,AMK)
       FORMA1 = FA1A1P(QQ)
       FORM1 = FORMA1*FORMKS

      ELSEIF (MNUM.EQ.4) THEN
C ------------ pi0 pi0 K-  (K*-pi0)
       FORMKS = BWIGM(S1,AMKST,GAMKST,AMPI,AMK)
       FORMK1 = FK1AB(QQ,3)
       FORM1 = FORMK1*FORMKS

      ELSEIF (MNUM.EQ.5) THEN
C ------------ K- pi- pi+ (rho0 K-)
       FORMK1 = FK1AB(QQ,4)
       FORMRO = FPIKM(SQRT(S1),AMPI,AMPI)
       FORM1 = FORMK1*FORMRO

      ELSEIF (MNUM.EQ.6) THEN
C ------------ pi- K0B pi0 (pi- K*0B)
       FORMK1 = FK1AB(QQ,1)
       FORMKS = BWIGM(S1,AMKST,GAMKST,AMK,AMPI)
       FORM1 = FORMK1*FORMKS

      ELSEIF (MNUM.EQ.7) THEN
C -------------- eta pi- pi0 final state
       FORM1=0.0
      ENDIF
      END
      FUNCTION FORM2(MNUM,QQ,S1,SDWA)
C     ==================================================================
C     formfactorfor F2 for 3 scalar final state
C     R. Fisher, J. Wess and F. Wagner Z. Phys C3 (1980) 313
C     H. Georgi, Weak interactions and modern particle theory,
C     The Benjamin/Cummings Pub. Co., Inc. 1984.
C     R. Decker, E. Mirkes, R. Sauer, Z. Was Karlsruhe preprint TTP92-25
C     and erratum !!!!!!
C     ==================================================================
C
      COMPLEX FORM2,WIGNER,WIGFOR,FPIKM,BWIGM
      COMMON / PARMAS / AMTAU,AMNUTA,AMEL,AMNUE,AMMU,AMNUMU
     *                 ,AMPIZ,AMPI,AMRO,GAMRO,AMA1,GAMA1
     *                 ,AMK,AMKZ,AMKST,GAMKST
C
      REAL*4            AMTAU,AMNUTA,AMEL,AMNUE,AMMU,AMNUMU
     *                 ,AMPIZ,AMPI,AMRO,GAMRO,AMA1,GAMA1
     *                 ,AMK,AMKZ,AMKST,GAMKST
      COMMON /IPChT/ IVER
      INTEGER        IVER

      COMPLEX FORMA1,FORMK1,FORMRO,FORMKS
      COMPLEX FA1A1P,FK1AB,F3PI,F3PI_RCHT

      IF     (MNUM.EQ.0) THEN
C ------------  3 pi hadronic state (a1)
C       FORMRO = FPIKM(SQRT(S1),AMPI,AMPI)
C       FORMRO = F3PI(2,QQ,S1,SDWA)
C       FORMA1 = FA1A1P(QQ)
C       FORM2 = FORMA1*FORMRO
C       FORM2 = F3PI(2,QQ,S1,SDWA)
      IF (IVER.EQ.0) THEN
       FORM2 = F3PI(2,QQ,S1,SDWA)
      ELSE
       FORM2 = F3PI_RCHT(2,QQ,S1,SDWA)
      ENDIF
      ELSEIF (MNUM.EQ.1) THEN
C ------------ K- pi- K+ (rho0 pi-)
       FORMRO = FPIKM(SQRT(S1),AMK,AMK)
       FORMA1 = FA1A1P(QQ)
       FORM2 = FORMA1*FORMRO

      ELSEIF (MNUM.EQ.2) THEN
C ------------ K0 pi- K0B (rho0 pi-)
       FORMRO = FPIKM(SQRT(S1),AMK,AMK)
       FORMA1 = FA1A1P(QQ)
       FORM2 = FORMA1*FORMRO

      ELSEIF (MNUM.EQ.3) THEN
C ------------ K- pi0 K0 (rho- pi0)
       FORMRO = FPIKM(SQRT(S1),AMK,AMK)
       FORMA1 = FA1A1P(QQ)
       FORM2 = FORMA1*FORMRO

      ELSEIF (MNUM.EQ.4) THEN
C ------------ pi0 pi0 K-  (K*-pi0)
       FORMKS = BWIGM(S1,AMKST,GAMKST,AMPI,AMK)
       FORMK1 = FK1AB(QQ,3)
       FORM2 = FORMK1*FORMKS

      ELSEIF (MNUM.EQ.5) THEN
C ------------ K- pi- pi+  (K*0B pi-)
       FORMKS = BWIGM(S1,AMKST,GAMKST,AMPI,AMK)
       FORMK1 = FK1AB(QQ,1)
       FORM2 = FORMK1*FORMKS
C
      ELSEIF (MNUM.EQ.6) THEN
C ------------ pi- K0B pi0 (rho- K0B)
       FORMRO = FPIKM(SQRT(S1),AMPI,AMPI)
       FORMK1 = FK1AB(QQ,2)
       FORM2 = FORMK1*FORMRO
C
      ELSEIF (MNUM.EQ.7) THEN
C -------------- eta pi- pi0 final state
       FORM2=0.0
      ENDIF
C
      END
      COMPLEX FUNCTION BWIGM(S,M,G,XM1,XM2)
C **********************************************************
C     P-WAVE BREIT-WIGNER  FOR RHO
C **********************************************************
      REAL S,M,G,XM1,XM2
      REAL PI,QS,QM,W,GS
      DATA INIT /0/
C ------------ PARAMETERS --------------------
      IF (INIT.EQ.0) THEN
      INIT=1
      PI=3.141592654
C -------  BREIT-WIGNER -----------------------
         ENDIF
       IF (S.GT.(XM1+XM2)**2) THEN
         QS=SQRT(ABS((S   -(XM1+XM2)**2)*(S   -(XM1-XM2)**2)))/SQRT(S)
         QM=SQRT(ABS((M**2-(XM1+XM2)**2)*(M**2-(XM1-XM2)**2)))/M
         W=SQRT(S)
         GS=G*(M/W)**2*(QS/QM)**3
       ELSE
         GS=0.0
       ENDIF
         BWIGM=M**2/CMPLX(M**2-S,-SQRT(S)*GS)
      RETURN
      END
      COMPLEX FUNCTION FPIKM(W,XM1,XM2)
C **********************************************************
C     PION FORM FACTOR
C **********************************************************
      COMPLEX BWIGM
      REAL ROM,ROG,ROM1,ROG1,BETA1,PI,PIM,S,W
      EXTERNAL BWIG
      DATA  INIT /0/
C
C ------------ PARAMETERS --------------------
      IF (INIT.EQ.0 ) THEN
      INIT=1
      PI=3.141592654
      PIM=.140
      ROM=0.773
      ROG=0.145
      ROM1=1.370
      ROG1=0.510
      BETA1=-0.145
      ENDIF
C -----------------------------------------------
      S=W**2
      FPIKM=(BWIGM(S,ROM,ROG,XM1,XM2)+BETA1*BWIGM(S,ROM1,ROG1,XM1,XM2))
     & /(1+BETA1)
      RETURN
      END
      COMPLEX FUNCTION FPIKMD(W,XM1,XM2)
C **********************************************************
C     PION FORM FACTOR
C **********************************************************
      COMPLEX BWIGM
      REAL ROM,ROG,ROM1,ROG1,PI,PIM,S,W
      EXTERNAL BWIG
      DATA  INIT /0/
C
C ------------ PARAMETERS --------------------
      IF (INIT.EQ.0 ) THEN
      INIT=1
      PI=3.141592654
      PIM=.140
      ROM=0.773
      ROG=0.145
      ROM1=1.500
      ROG1=0.220
      ROM2=1.750
      ROG2=0.120
      BETA=6.5
      DELTA=-26.0
      ENDIF
C -----------------------------------------------
      S=W**2
      FPIKMD=(DELTA*BWIGM(S,ROM,ROG,XM1,XM2)
     $      +BETA*BWIGM(S,ROM1,ROG1,XM1,XM2)
     $      +     BWIGM(S,ROM2,ROG2,XM1,XM2))
     & /(1+BETA+DELTA)
      RETURN
      END
 
      FUNCTION FORM3(MNUM,QQ,S1,SDWA)
C     ==================================================================
C     formfactorfor F3 for 3 scalar final state
C     R. Fisher, J. Wess and F. Wagner Z. Phys C3 (1980) 313
C     H. Georgi, Weak interactions and modern particle theory,
C     The Benjamin/Cummings Pub. Co., Inc. 1984.
C     R. Decker, E. Mirkes, R. Sauer, Z. Was Karlsruhe preprint TTP92-25
C     and erratum !!!!!!
C     ==================================================================
C
      COMMON / PARMAS / AMTAU,AMNUTA,AMEL,AMNUE,AMMU,AMNUMU
     *                 ,AMPIZ,AMPI,AMRO,GAMRO,AMA1,GAMA1
     *                 ,AMK,AMKZ,AMKST,GAMKST
C
      REAL*4            AMTAU,AMNUTA,AMEL,AMNUE,AMMU,AMNUMU
     *                 ,AMPIZ,AMPI,AMRO,GAMRO,AMA1,GAMA1
     *                 ,AMK,AMKZ,AMKST,GAMKST
      COMMON /IPChT/ IVER
      INTEGER        IVER
      COMPLEX FORM3,BWIGM
      COMPLEX FORMA1,FORMK1,FORMRO,FORMKS
      COMPLEX FA1A1P,FK1AB,F3PI,F3PI_RCHT
C
      IF (MNUM.EQ.0) THEN
C ------------  3 pi hadronic state (a1)
C       FORMRO = FPIKM(SQRT(S1),AMPI,AMPI)
C       FORMRO = F3PI(3,QQ,S1,SDWA)
C       FORMA1 = FA1A1P(QQ)
C       FORM3 = FORMA1*FORMRO
      IF (IVER.EQ.0) THEN
        FORM3 = F3PI(3,QQ,S1,SDWA)
      ELSE
        FORM3 = (0.,0.)
      ENDIF

      ELSEIF (MNUM.EQ.3) THEN
C ------------ K- pi0 K0  (K*- K0)
       FORMKS = BWIGM(S1,AMKST,GAMKST,AMPIZ,AMK)
       FORMA1 = FA1A1P(QQ)
       FORM3 = FORMA1*FORMKS

      ELSEIF (MNUM.EQ.6) THEN
C ------------ pi- K0B pi0 (K*- pi0)
       FORMKS = BWIGM(S1,AMKST,GAMKST,AMK,AMPI)
       FORMK1 = FK1AB(QQ,3)
       FORM3 = FORMK1*FORMKS

      ELSE
       FORM3=CMPLX(0.,0.)
      ENDIF
      END
      FUNCTION FORM4(MNUM,QQ,S1,S2,S3)
C     ==================================================================
C     formfactorfor F4 for 3 scalar final state
C     R. Decker, in preparation
C     R. Decker, E. Mirkes, R. Sauer, Z. Was Karlsruhe preprint TTP92-25
C     and erratum !!!!!!
C     ==================================================================
C
      COMMON / PARMAS / AMTAU,AMNUTA,AMEL,AMNUE,AMMU,AMNUMU
     *                 ,AMPIZ,AMPI,AMRO,GAMRO,AMA1,GAMA1
     *                 ,AMK,AMKZ,AMKST,GAMKST
C
      REAL*4            AMTAU,AMNUTA,AMEL,AMNUE,AMMU,AMNUMU
     *                 ,AMPIZ,AMPI,AMRO,GAMRO,AMA1,GAMA1
     *                 ,AMK,AMKZ,AMKST,GAMKST
      COMMON /IPChT/ IVER
      INTEGER        IVER

      COMPLEX FORM4,WIGNER,FPIKM,F3PI_RCHT
      REAL*4 M
C ---- this formfactor is switched off for cleo version
       FORM4=CMPLX(0.0,0.0)
      IF (MNUM.EQ.0) THEN
C ------------  3 pi hadronic state (a1)
C       FORMRO = FPIKM(SQRT(S1),AMPI,AMPI)
C       FORMRO = F3PI(3,QQ,S1,SDWA)
C       FORMA1 = FA1A1P(QQ)
C       FORM3 = FORMA1*FORMRO
C        FORM4 = CMPLX(-1., 0.) 
       IF (IVER.EQ.1) FORM4 = F3PI_RCHT(4,QQ,S1,S2)* CMPLX(0., 1.)
      ENDIF


      END
      FUNCTION FORM5(MNUM,QQ,S1,S2)
C     ==================================================================
C     formfactorfor F5 for 3 scalar final state
C     G. Kramer, W. Palmer, S. Pinsky, Phys. Rev. D30 (1984) 89.
C     G. Kramer, W. Palmer             Z. Phys. C25 (1984) 195.
C     R. Decker, E. Mirkes, R. Sauer, Z. Was Karlsruhe preprint TTP92-25
C     and erratum !!!!!!
C     ==================================================================
C
      COMMON / PARMAS / AMTAU,AMNUTA,AMEL,AMNUE,AMMU,AMNUMU
     *                 ,AMPIZ,AMPI,AMRO,GAMRO,AMA1,GAMA1
     *                 ,AMK,AMKZ,AMKST,GAMKST
C
      REAL*4            AMTAU,AMNUTA,AMEL,AMNUE,AMMU,AMNUMU
     *                 ,AMPIZ,AMPI,AMRO,GAMRO,AMA1,GAMA1
     *                 ,AMK,AMKZ,AMKST,GAMKST
      COMPLEX FORM5,WIGNER,FPIKM,FPIKMD,BWIGM
      IF     (MNUM.EQ.0) THEN
C ------------  3 pi hadronic state (a1)
        FORM5=0.0
      ELSEIF (MNUM.EQ.1) THEN
C ------------ K- pi- K+
         ELPHA=-0.2
         FORM5=FPIKMD(SQRT(QQ),AMPI,AMPI)/(1+ELPHA)
     $        *(       FPIKM(SQRT(S2),AMPI,AMPI)
     $          +ELPHA*BWIGM(S1,AMKST,GAMKST,AMPI,AMK))
      ELSEIF (MNUM.EQ.2) THEN
C ------------ K0 pi- K0B
         ELPHA=-0.2
         FORM5=FPIKMD(SQRT(QQ),AMPI,AMPI)/(1+ELPHA)
     $        *(       FPIKM(SQRT(S2),AMPI,AMPI)
     $          +ELPHA*BWIGM(S1,AMKST,GAMKST,AMPI,AMK))
      ELSEIF (MNUM.EQ.3) THEN
C ------------ K- K0 pi0
        FORM5=0.0
      ELSEIF (MNUM.EQ.4) THEN
C ------------ pi0 pi0 K-
        FORM5=0.0
      ELSEIF (MNUM.EQ.5) THEN
C ------------ K- pi- pi+
        ELPHA=-0.2
        FORM5=BWIGM(QQ,AMKST,GAMKST,AMPI,AMK)/(1+ELPHA)
     $       *(       FPIKM(SQRT(S1),AMPI,AMPI)
     $         +ELPHA*BWIGM(S2,AMKST,GAMKST,AMPI,AMK))
      ELSEIF (MNUM.EQ.6) THEN
C ------------ pi- K0B pi0
        ELPHA=-0.2
        FORM5=BWIGM(QQ,AMKST,GAMKST,AMPI,AMKZ)/(1+ELPHA)
     $       *(       FPIKM(SQRT(S2),AMPI,AMPI)
     $         +ELPHA*BWIGM(S1,AMKST,GAMKST,AMPI,AMK))
      ELSEIF (MNUM.EQ.7) THEN
C -------------- eta pi- pi0 final state
       FORM5=FPIKMD(SQRT(QQ),AMPI,AMPI)*FPIKM(SQRT(S1),AMPI,AMPI)
      ENDIF
C
      END
      SUBROUTINE CURRX(MNUM,PIM1,PIM2,PIM3,PIM4,HADCUR)
C     ==================================================================
C     hadronic current for 4 pi final state
C     R. Fisher, J. Wess and F. Wagner Z. Phys C3 (1980) 313
C     R. Decker Z. Phys C36 (1987) 487.
C     M. Gell-Mann, D. Sharp, W. Wagner Phys. Rev. Lett 8 (1962) 261.
C     ==================================================================
 
      COMMON / PARMAS / AMTAU,AMNUTA,AMEL,AMNUE,AMMU,AMNUMU
     *                 ,AMPIZ,AMPI,AMRO,GAMRO,AMA1,GAMA1
     *                 ,AMK,AMKZ,AMKST,GAMKST
C
      REAL*4            AMTAU,AMNUTA,AMEL,AMNUE,AMMU,AMNUMU
     *                 ,AMPIZ,AMPI,AMRO,GAMRO,AMA1,GAMA1
     *                 ,AMK,AMKZ,AMKST,GAMKST
      COMMON / DECPAR / GFERMI,GV,GA,CCABIB,SCABIB,GAMEL
      REAL*4            GFERMI,GV,GA,CCABIB,SCABIB,GAMEL
C ARBITRARY FIXING OF THE FOUR PI X-SECTION NORMALIZATION
      COMMON /ARBIT/ ARFLAT,AROMEG
      REAL  PIM1(4),PIM2(4),PIM3(4),PIM4(4),PAA(4)
      COMPLEX HADCUR(4),FORM1,FORM2,FORM3,FPIKM
      COMPLEX BWIGN
      REAL PA(4),PB(4)
      REAL AA(4,4),PP(4,4)
      DATA PI /3.141592653589793238462643/
      DATA  FPI /93.3E-3/
      BWIGN(A,XM,XG)=1.0/CMPLX(A-XM**2,XM*XG)
C
C --- masses and constants
      G1=12.924
      G2=1475.98
      G =G1*G2
      ELPHA=-.1
      AMROP=1.7
      GAMROP=0.26
      AMOM=.782
      GAMOM=0.0085
      ARFLAT=1.0
      AROMEG=1.0
C
      FRO=0.266*AMRO**2
      COEF1=2.0*SQRT(3.0)/FPI**2*ARFLAT
      COEF2=FRO*G*AROMEG
C --- initialization of four vectors
      DO 7 K=1,4
      DO 8 L=1,4
 8    AA(K,L)=0.0
      HADCUR(K)=CMPLX(0.0)
      PAA(K)=PIM1(K)+PIM2(K)+PIM3(K)+PIM4(K)
      PP(1,K)=PIM1(K)
      PP(2,K)=PIM2(K)
      PP(3,K)=PIM3(K)
 7    PP(4,K)=PIM4(K)
C
      IF (MNUM.EQ.1) THEN
C ===================================================================
C pi- pi- p0 pi+ case                                            ====
C ===================================================================
       QQ=PAA(4)**2-PAA(3)**2-PAA(2)**2-PAA(1)**2
C --- loop over thre contribution of the non-omega current
       DO 201 K=1,3
        SK=(PP(K,4)+PIM4(4))**2-(PP(K,3)+PIM4(3))**2
     $    -(PP(K,2)+PIM4(2))**2-(PP(K,1)+PIM4(1))**2
C -- definition of AA matrix
C -- cronecker delta
        DO 202 I=1,4
         DO 203 J=1,4
 203     AA(I,J)=0.0
 202    AA(I,I)=1.0
C ... and the rest ...
        DO 204 L=1,3
         IF (L.NE.K) THEN
          DENOM=(PAA(4)-PP(L,4))**2-(PAA(3)-PP(L,3))**2
     $         -(PAA(2)-PP(L,2))**2-(PAA(1)-PP(L,1))**2
          DO 205 I=1,4
          DO 205 J=1,4
                      SIG= 1.0
           IF(J.NE.4) SIG=-SIG
           AA(I,J)=AA(I,J)
     $            -SIG*(PAA(I)-2.0*PP(L,I))*(PAA(J)-PP(L,J))/DENOM
 205      CONTINUE
         ENDIF
 204    CONTINUE
C --- lets add something to HADCURR
       FORM1= FPIKM(SQRT(SK),AMPI,AMPI) *FPIKM(SQRT(QQ),AMPI,AMPI)
C       FORM1= FPIKM(SQRT(SK),AMPI,AMPI) *FPIKMD(SQRT(QQ),AMPI,AMPI)
CCCCCCCCCCCCCCCCC       FORM1=WIGFOR(SK,AMRO,GAMRO)      (tests)
C
       FIX=1.0
       IF (K.EQ.3) FIX=-2.0
       DO 206 I=1,4
       DO 206 J=1,4
        HADCUR(I)=
     $  HADCUR(I)+CMPLX(FIX*COEF1)*FORM1*AA(I,J)*(PP(K,J)-PP(4,J))
 206   CONTINUE
C --- end of the non omega current (3 possibilities)
 201   CONTINUE
C
C
C --- there are two possibilities for omega current
C --- PA PB are corresponding first and second pi-s
       DO 301 KK=1,2
        DO 302 I=1,4
         PA(I)=PP(KK,I)
         PB(I)=PP(3-KK,I)
 302    CONTINUE
C --- lorentz invariants
         QQA=0.0
         SS23=0.0
         SS24=0.0
         SS34=0.0
         QP1P2=0.0
         QP1P3=0.0
         QP1P4=0.0
         P1P2 =0.0
         P1P3 =0.0
         P1P4 =0.0
        DO 303 K=1,4
                     SIGN=-1.0
         IF (K.EQ.4) SIGN= 1.0
         QQA=QQA+SIGN*(PAA(K)-PA(K))**2
         SS23=SS23+SIGN*(PB(K)  +PIM3(K))**2
         SS24=SS24+SIGN*(PB(K)  +PIM4(K))**2
         SS34=SS34+SIGN*(PIM3(K)+PIM4(K))**2
         QP1P2=QP1P2+SIGN*(PAA(K)-PA(K))*PB(K)
         QP1P3=QP1P3+SIGN*(PAA(K)-PA(K))*PIM3(K)
         QP1P4=QP1P4+SIGN*(PAA(K)-PA(K))*PIM4(K)
         P1P2=P1P2+SIGN*PA(K)*PB(K)
         P1P3=P1P3+SIGN*PA(K)*PIM3(K)
         P1P4=P1P4+SIGN*PA(K)*PIM4(K)
 303    CONTINUE
C
        FORM2=COEF2*(BWIGN(QQ,AMRO,GAMRO)+ELPHA*BWIGN(QQ,AMROP,GAMROP))
C        FORM3=BWIGN(QQA,AMOM,GAMOM)*(BWIGN(SS23,AMRO,GAMRO)+
C     $        BWIGN(SS24,AMRO,GAMRO)+BWIGN(SS34,AMRO,GAMRO))
        FORM3=BWIGN(QQA,AMOM,GAMOM)
C
        DO 304 K=1,4
         HADCUR(K)=HADCUR(K)+FORM2*FORM3*(
     $             PB  (K)*(QP1P3*P1P4-QP1P4*P1P3)
     $            +PIM3(K)*(QP1P4*P1P2-QP1P2*P1P4)
     $            +PIM4(K)*(QP1P2*P1P3-QP1P3*P1P2) )
 304    CONTINUE
 301   CONTINUE
C
      ELSE
C ===================================================================
C pi0 pi0 p0 pi- case                                            ====
C ===================================================================
       QQ=PAA(4)**2-PAA(3)**2-PAA(2)**2-PAA(1)**2
       DO 101 K=1,3
C --- loop over thre contribution of the non-omega current
        SK=(PP(K,4)+PIM4(4))**2-(PP(K,3)+PIM4(3))**2
     $    -(PP(K,2)+PIM4(2))**2-(PP(K,1)+PIM4(1))**2
C -- definition of AA matrix
C -- cronecker delta
        DO 102 I=1,4
         DO 103 J=1,4
 103     AA(I,J)=0.0
 102    AA(I,I)=1.0
C
C ... and the rest ...
        DO 104 L=1,3
         IF (L.NE.K) THEN
          DENOM=(PAA(4)-PP(L,4))**2-(PAA(3)-PP(L,3))**2
     $         -(PAA(2)-PP(L,2))**2-(PAA(1)-PP(L,1))**2
          DO 105 I=1,4
          DO 105 J=1,4
                      SIG=1.0
           IF(J.NE.4) SIG=-SIG
           AA(I,J)=AA(I,J)
     $            -SIG*(PAA(I)-2.0*PP(L,I))*(PAA(J)-PP(L,J))/DENOM
 105      CONTINUE
         ENDIF
 104    CONTINUE
C --- lets add something to HADCURR
       FORM1= FPIKM(SQRT(SK),AMPI,AMPI) *FPIKM(SQRT(QQ),AMPI,AMPI)
C       FORM1= FPIKM(SQRT(SK),AMPI,AMPI) *FPIKMD(SQRT(QQ),AMPI,AMPI)
CCCCCCCCCCCCC       FORM1=WIGFOR(SK,AMRO,GAMRO)        (tests)
        DO 106 I=1,4
        DO 106 J=1,4
         HADCUR(I)=
     $   HADCUR(I)+CMPLX(COEF1)*FORM1*AA(I,J)*(PP(K,J)-PP(4,J))
 106    CONTINUE
C --- end of the non omega current (3 possibilities)
 101   CONTINUE
      ENDIF
      END
      FUNCTION WIGFOR(S,XM,XGAM)
      COMPLEX WIGFOR,WIGNOR
      WIGNOR=CMPLX(-XM**2,XM*XGAM)
      WIGFOR=WIGNOR/CMPLX(S-XM**2,XM*XGAM)
      END
      SUBROUTINE CURINF
C HERE the form factors of M. Finkemeier et al. start
C it ends with the string:  M. Finkemeier et al. END
      COMMON /INOUT/ INUT, IOUT
      WRITE (UNIT = IOUT,FMT = 99)
      WRITE (UNIT = IOUT,FMT = 98)
c                    print *, 'here is curinf'
 99   FORMAT(
     . /,   ' *************************************************** ',
     . /,   '   YOU ARE USING THE 4 PION DECAY MODE FORM FACTORS    ',
     . /,   '   WHICH HAVE BEEN DESCRIBED IN:',
     . /,   ' R. DECKER, M. FINKEMEIER, P. HEILIGER AND H.H. JONSSON',
     . /,   '   "TAU DECAYS INTO FOUR PIONS" ',
     . /,   '   UNIVERSITAET KARLSRUHE PREPRINT TTP 94-13 (1994);',
     . /,   '                    LNF-94/066(IR); HEP-PH/9410260  ',
     . /,   '  ',
     . /,   ' PLEASE NOTE THAT THIS ROUTINE IS USING PARAMETERS',
     . /,   ' RELATED TO THE 3 PION DECAY MODE (A1 MODE), SUCH AS',
     . /,   ' THE A1 MASS AND WIDTH (TAKEN FROM THE COMMON /PARMAS/)',
     . /,   ' AND THE 2 PION VECTOR RESONANCE FORM FACTOR (BY USING',
     . /,   ' THE ROUTINE FPIKM)'                                   ,
     . /,   ' THUS IF YOU DECIDE TO CHANGE ANY OF THESE, YOU WILL'  ,
     . /,   ' HAVE TO REFIT THE 4 PION PARAMETERS IN THE COMMON'    )
   98   FORMAT(
     .      ' BLOCK /TAU4PI/, OR YOU MIGHT GET A BAD DISCRIPTION'   ,
     . /,   ' OF TAU -> 4 PIONS'       ,
     . /,   ' for these formfactors set in routine CHOICE for',
     . /,  ' mnum.eq.102 -- AMRX=1.42 and GAMRX=.21',
     . /,  ' mnum.eq.101 -- AMRX=1.3 and GAMRX=.46 PROB1,PROB2=0.2',
     . /,  ' to optimize phase space parametrization',
     . /,   ' *************************************************** ',
     . /,   ' coded by M. Finkemeier and P. Heiliger, 29. sept. 1994',
     . /,   ' incorporated to TAUOLA by Z. Was      17. jan. 1995',
c     . /,   ' fitted on (day/month/year) by ...  ',
c     . /,   ' to .... data ',
     . /,   ' changed by: Z. Was on 17.01.95',
     . /,   ' changes by: M. Finkemeier on 30.01.95' )
      END
C
      SUBROUTINE CURINI
      COMMON /TAU4PI/ GOMEGA,GAMMA1,GAMMA2,ROM1,ROG1,BETA1,
     .                ROM2,ROG2,BETA2
      REAL*4          GOMEGA,GAMMA1,GAMMA2,ROM1,ROG1,BETA1,
     .                ROM2,ROG2,BETA2
      GOMEGA = 1.4
      GAMMA1 = 0.38
      GAMMA2 = 0.38
      ROM1   = 1.35
      ROG1   = 0.3
      BETA1  = 0.08
      ROM2   = 1.70
      ROG2   = 0.235
      BETA2  = -0.0075				
      END						
      COMPLEX FUNCTION BWIGA1(QA)
C     ================================================================
C     breit-wigner enhancement of a1
C     ================================================================
      COMPLEX WIGNER
      COMMON / PARMAS/ AMTAU,AMNUTA,AMEL,AMNUE,AMMU,AMNUMU,
     %                 AMPIZ,AMPI,AMRO,GAMRO,AMA1,GAMA1,
     %                 AMK,AMKZ,AMKST,GAMKST
 
C
      REAL*4           AMTAU,AMNUTA,AMEL,AMNUE,AMMU,AMNUMU,
     %                 AMPIZ,AMPI,AMRO,GAMRO,AMA1,GAMA1,
     %                 AMK,AMKZ,AMKST,GAMKST
 
      WIGNER(A,B,C)=CMPLX(1.0,0.0)/CMPLX(A-B**2,B*C)
      GAMAX=GAMA1*GFUN(QA)/GFUN(AMA1**2)
      BWIGA1=-AMA1**2*WIGNER(QA,AMA1,GAMAX)
      RETURN
      END
      COMPLEX FUNCTION BWIGEPS(QEPS)
C     =============================================================
C     breit-wigner enhancement of epsilon
C     =============================================================
      REAL QEPS,MEPS,GEPS
      MEPS=1.300
      GEPS=.600
      BWIGEPS=CMPLX(MEPS**2,-MEPS*GEPS)/
     %        CMPLX(MEPS**2-QEPS,-MEPS*GEPS)
      RETURN
      END
      COMPLEX FUNCTION FRHO4(W,XM1,XM2)
C     ===========================================================
C     rho-type resonance factor with higher radials, to be used
C     by CURR for the four pion mode
C     ===========================================================
      COMPLEX BWIGM
      COMMON /TAU4PI/ GOMEGA,GAMMA1,GAMMA2,ROM1,ROG1,BETA1,
     .                ROM2,ROG2,BETA2
      REAL*4          GOMEGA,GAMMA1,GAMMA2,ROM1,ROG1,BETA1,
     .                ROM2,ROG2,BETA2
      REAL ROM,ROG,PI,PIM,S,W
      EXTERNAL BWIG
      DATA  INIT /0/
C
C ------------ PARAMETERS --------------------
      IF (INIT.EQ.0 ) THEN
      INIT=1
      PI=3.141592654
      PIM=.140
      ROM=0.773
      ROG=0.145
      ENDIF
C -----------------------------------------------
      S=W**2
c	       print *,'rom2,rog2 =',rom2,rog2
      FRHO4=(BWIGM(S,ROM,ROG,XM1,XM2)+BETA1*BWIGM(S,ROM1,ROG1,XM1,XM2)
     & +BETA2*BWIGM(S,ROM2,ROG2,XM1,XM2))
     & /(1+BETA1+BETA2)
      RETURN
      END
      SUBROUTINE CURR(MNUM,PIM1,PIM2,PIM3,PIM4,HADCUR)
C     ==================================================================
C     Hadronic current for 4 pi final state, according to:
C     R. Decker, M. Finkemeier, P. Heiliger, H.H.Jonsson, TTP94-13
C
C     See also:
C     R. Fisher, J. Wess and F. Wagner Z. Phys C3 (1980) 313
C     R. Decker Z. Phys C36 (1987) 487.
C     M. Gell-Mann, D. Sharp, W. Wagner Phys. Rev. Lett 8 (1962) 261.
C     ==================================================================
 
      COMMON /TAU4PI/ GOMEGA,GAMMA1,GAMMA2,ROM1,ROG1,BETA1,
     .                ROM2,ROG2,BETA2
      REAL*4          GOMEGA,GAMMA1,GAMMA2,ROM1,ROG1,BETA1,
     .                ROM2,ROG2,BETA2
      COMMON / PARMAS / AMTAU,AMNUTA,AMEL,AMNUE,AMMU,AMNUMU
     *                 ,AMPIZ,AMPI,AMRO,GAMRO,AMA1,GAMA1
     *                 ,AMK,AMKZ,AMKST,GAMKST
C
      REAL*4            AMTAU,AMNUTA,AMEL,AMNUE,AMMU,AMNUMU
     *                 ,AMPIZ,AMPI,AMRO,GAMRO,AMA1,GAMA1
     *                 ,AMK,AMKZ,AMKST,GAMKST
      COMMON / DECPAR / GFERMI,GV,GA,CCABIB,SCABIB,GAMEL
      REAL*4            GFERMI,GV,GA,CCABIB,SCABIB,GAMEL
      REAL  PIM1(4),PIM2(4),PIM3(4),PIM4(4),PAA(4)
      COMPLEX HADCUR(4),FORM1,FORM2,FORM3,FPIKM
      COMPLEX BWIGN,FRHO4
      COMPLEX BWIGEPS,BWIGA1
      COMPLEX HCOMP1(4),HCOMP2(4),HCOMP3(4),HCOMP4(4)
 
      COMPLEX T243,T213,T143,T123,T341,T342
      COMPLEX T124,T134,T214,T234,T314,T324
      COMPLEX S2413,S2314,S1423,S1324,S34
      COMPLEX S2431,S3421
      COMPLEX BRACK1,BRACK2,BRACK3,BRACK4A,BRACK4B,BRACK4
 
      REAL QMP1,QMP2,QMP3,QMP4
      REAL PS43,PS41,PS42,PS34,PS14,PS13,PS24,PS23
      REAL PS21,PS31
 
      REAL PD243,PD241,PD213,PD143,PD142
      REAL PD123,PD341,PD342,PD413,PD423
      REAL PD124,PD134,PD214,PD234,PD314,PD324
      REAL QP1,QP2,QP3,QP4
 
      REAL PA(4),PB(4)
      REAL AA(4,4),PP(4,4)
      DATA PI /3.141592653589793238462643/
      DATA  FPI /93.3E-3/
      DATA INIT /0/	
      BWIGN(A,XM,XG)=1.0/CMPLX(A-XM**2,XM*XG)
C
      IF (INIT.EQ.0) THEN
	CALL CURINI
	CALL CURINF
	INIT = 1
      ENDIF	        	
C
C --- MASSES AND CONSTANTS
      G1=12.924
      G2=1475.98 * GOMEGA
      G =G1*G2
      ELPHA=-.1
      AMROP=1.7
      GAMROP=0.26
      AMOM=.782
      GAMOM=0.0085
      ARFLAT=1.0
      AROMEG=1.0
C
      FRO=0.266*AMRO**2
      COEF1=2.0*SQRT(3.0)/FPI**2*ARFLAT
      COEF2=FRO*G*AROMEG
C --- INITIALIZATION OF FOUR VECTORS
      DO 7 K=1,4
      DO 8 L=1,4
 8    AA(K,L)=0.0
      HADCUR(K)=CMPLX(0.0)
      PAA(K)=PIM1(K)+PIM2(K)+PIM3(K)+PIM4(K)
      PP(1,K)=PIM1(K)
      PP(2,K)=PIM2(K)
      PP(3,K)=PIM3(K)
 7    PP(4,K)=PIM4(K)
C
      IF (MNUM.EQ.1) THEN
C ===================================================================
C PI- PI- P0 PI+ CASE                                            ====
C ===================================================================
       QQ=PAA(4)**2-PAA(3)**2-PAA(2)**2-PAA(1)**2
 
C FIRST DEFINITION OF SCALAR PRODUCTS OF MOMENTUM VECTORS
 
C DEFINE (Q-PI)**2 AS QPI:
 
      QMP1=(PIM2(4)+PIM3(4)+PIM4(4))**2-(PIM2(3)+PIM3(3)+PIM4(3))**2
     %   -(PIM2(2)+PIM3(2)+PIM4(2))**2-(PIM2(1)+PIM3(1)+PIM4(1))**2
 
      QMP2=(PIM1(4)+PIM3(4)+PIM4(4))**2-(PIM1(3)+PIM3(3)+PIM4(3))**2
     %   -(PIM1(2)+PIM3(2)+PIM4(2))**2-(PIM1(1)+PIM3(1)+PIM4(1))**2
 
      QMP3=(PIM1(4)+PIM2(4)+PIM4(4))**2-(PIM1(3)+PIM2(3)+PIM4(3))**2
     %   -(PIM1(2)+PIM2(2)+PIM4(2))**2-(PIM1(1)+PIM2(1)+PIM4(1))**2
 
      QMP4=(PIM1(4)+PIM2(4)+PIM3(4))**2-(PIM1(3)+PIM2(3)+PIM3(3))**2
     %   -(PIM1(2)+PIM2(2)+PIM3(2))**2-(PIM1(1)+PIM2(1)+PIM3(1))**2
 
 
C DEFINE (PI+PK)**2 AS PSIK:
 
      PS43=(PIM4(4)+PIM3(4))**2-(PIM4(3)+PIM3(3))**2
     %    -(PIM4(2)+PIM3(2))**2-(PIM4(1)+PIM3(1))**2
 
      PS41=(PIM4(4)+PIM1(4))**2-(PIM4(3)+PIM1(3))**2
     %    -(PIM4(2)+PIM1(2))**2-(PIM4(1)+PIM1(1))**2
 
      PS42=(PIM4(4)+PIM2(4))**2-(PIM4(3)+PIM2(3))**2
     %    -(PIM4(2)+PIM2(2))**2-(PIM4(1)+PIM2(1))**2
 
      PS34=PS43
 
      PS14=PS41
 
      PS13=(PIM1(4)+PIM3(4))**2-(PIM1(3)+PIM3(3))**2
     %    -(PIM1(2)+PIM3(2))**2-(PIM1(1)+PIM3(1))**2
 
      PS24=PS42
 
      PS23=(PIM2(4)+PIM3(4))**2-(PIM2(3)+PIM3(3))**2
     %    -(PIM2(2)+PIM3(2))**2-(PIM2(1)+PIM3(1))**2
 
      PD243=PIM2(4)*(PIM4(4)-PIM3(4))-PIM2(3)*(PIM4(3)-PIM3(3))
     %     -PIM2(2)*(PIM4(2)-PIM3(2))-PIM2(1)*(PIM4(1)-PIM3(1))
 
      PD241=PIM2(4)*(PIM4(4)-PIM1(4))-PIM2(3)*(PIM4(3)-PIM1(3))
     %     -PIM2(2)*(PIM4(2)-PIM1(2))-PIM2(1)*(PIM4(1)-PIM1(1))
 
      PD213=PIM2(4)*(PIM1(4)-PIM3(4))-PIM2(3)*(PIM1(3)-PIM3(3))
     %     -PIM2(2)*(PIM1(2)-PIM3(2))-PIM2(1)*(PIM1(1)-PIM3(1))
 
      PD143=PIM1(4)*(PIM4(4)-PIM3(4))-PIM1(3)*(PIM4(3)-PIM3(3))
     %     -PIM1(2)*(PIM4(2)-PIM3(2))-PIM1(1)*(PIM4(1)-PIM3(1))
 
      PD142=PIM1(4)*(PIM4(4)-PIM2(4))-PIM1(3)*(PIM4(3)-PIM2(3))
     %     -PIM1(2)*(PIM4(2)-PIM2(2))-PIM1(1)*(PIM4(1)-PIM2(1))
 
      PD123=PIM1(4)*(PIM2(4)-PIM3(4))-PIM1(3)*(PIM2(3)-PIM3(3))
     %     -PIM1(2)*(PIM2(2)-PIM3(2))-PIM1(1)*(PIM2(1)-PIM3(1))
 
      PD341=PIM3(4)*(PIM4(4)-PIM1(4))-PIM3(3)*(PIM4(3)-PIM1(3))
     %     -PIM3(2)*(PIM4(2)-PIM1(2))-PIM3(1)*(PIM4(1)-PIM1(1))
 
      PD342=PIM3(4)*(PIM4(4)-PIM2(4))-PIM3(3)*(PIM4(3)-PIM2(3))
     %     -PIM3(2)*(PIM4(2)-PIM2(2))-PIM3(1)*(PIM4(1)-PIM2(1))
 
      PD413=PIM4(4)*(PIM1(4)-PIM3(4))-PIM4(3)*(PIM1(3)-PIM3(3))
     %     -PIM4(2)*(PIM1(2)-PIM3(2))-PIM4(1)*(PIM1(1)-PIM3(1))
 
      PD423=PIM4(4)*(PIM2(4)-PIM3(4))-PIM4(3)*(PIM2(3)-PIM3(3))
     %     -PIM4(2)*(PIM2(2)-PIM3(2))-PIM4(1)*(PIM2(1)-PIM3(1))
 
C DEFINE Q*PI = QPI:
 
      QP1=PIM1(4)*(PIM1(4)+PIM2(4)+PIM3(4)+PIM4(4))
     %   -PIM1(3)*(PIM1(3)+PIM2(3)+PIM3(3)+PIM4(3))
     %   -PIM1(2)*(PIM1(2)+PIM2(2)+PIM3(2)+PIM4(2))
     %   -PIM1(1)*(PIM1(1)+PIM2(1)+PIM3(1)+PIM4(1))
 
      QP2=PIM2(4)*(PIM1(4)+PIM2(4)+PIM3(4)+PIM4(4))
     %   -PIM2(3)*(PIM1(3)+PIM2(3)+PIM3(3)+PIM4(3))
     %   -PIM2(2)*(PIM1(2)+PIM2(2)+PIM3(2)+PIM4(2))
     %   -PIM2(1)*(PIM1(1)+PIM2(1)+PIM3(1)+PIM4(1))
 
      QP3=PIM3(4)*(PIM1(4)+PIM2(4)+PIM3(4)+PIM4(4))
     %   -PIM3(3)*(PIM1(3)+PIM2(3)+PIM3(3)+PIM4(3))
     %   -PIM3(2)*(PIM1(2)+PIM2(2)+PIM3(2)+PIM4(2))
     %   -PIM3(1)*(PIM1(1)+PIM2(1)+PIM3(1)+PIM4(1))
 
      QP4=PIM4(4)*(PIM1(4)+PIM2(4)+PIM3(4)+PIM4(4))
     %   -PIM4(3)*(PIM1(3)+PIM2(3)+PIM3(3)+PIM4(3))
     %   -PIM4(2)*(PIM1(2)+PIM2(2)+PIM3(2)+PIM4(2))
     %   -PIM4(1)*(PIM1(1)+PIM2(1)+PIM3(1)+PIM4(1))
 
 
 
C DEFINE T(PI;PJ,PK)= TIJK:
 
      T243=BWIGA1(QMP2)*FPIKM(SQRT(PS43),AMPI,AMPI)*GAMMA1
      T213=BWIGA1(QMP2)*FPIKM(SQRT(PS13),AMPI,AMPI)*GAMMA1
      T143=BWIGA1(QMP1)*FPIKM(SQRT(PS43),AMPI,AMPI)*GAMMA1
      T123=BWIGA1(QMP1)*FPIKM(SQRT(PS23),AMPI,AMPI)*GAMMA1
      T341=BWIGA1(QMP3)*FPIKM(SQRT(PS41),AMPI,AMPI)*GAMMA1
      T342=BWIGA1(QMP3)*FPIKM(SQRT(PS42),AMPI,AMPI)*GAMMA1
 
C DEFINE S(I,J;K,L)= SIJKL:
 
      S2413=FRHO4(SQRT(PS24),AMPI,AMPI)*GAMMA2
      S2314=FRHO4(SQRT(PS23),AMPI,AMPI)*BWIGEPS(PS14)*GAMMA2
      S1423=FRHO4(SQRT(PS14),AMPI,AMPI)*GAMMA2
      S1324=FRHO4(SQRT(PS13),AMPI,AMPI)*BWIGEPS(PS24)*GAMMA2
      S34=FRHO4(SQRT(PS34),AMPI,AMPI)*GAMMA2
 
C DEFINITION OF AMPLITUDE, FIRST THE [] BRACKETS:
 
      BRACK1=2.*T143+2.*T243+T123+T213
     %    +T341*(PD241/QMP3-1.)+T342*(PD142/QMP3-1.)
     %    +3./4.*(S1423+S2413-S2314-S1324)-3.*S34
 
      BRACK2=2.*T143*PD243/QMP1+3.*T213
     %    +T123*(2.*PD423/QMP1+1.)+T341*(PD241/QMP3+3.)
     %    +T342*(PD142/QMP3+1.)
     %    -3./4.*(S2314+3.*S1324+3.*S1423+S2413)
 
      BRACK3=2.*T243*PD143/QMP2+3.*T123
     %    +T213*(2.*PD413/QMP2+1.)+T341*(PD241/QMP3+1.)
     %    +T342*(PD142/QMP3+3.)
     %    -3./4.*(3.*S2314+S1324+S1423+3.*S2413)
 
      BRACK4A=2.*T143*(PD243/QQ*(QP1/QMP1+1.)+PD143/QQ)
     %     +2.*T243*(PD143/QQ*(QP2/QMP2+1.)+PD243/QQ)
     %     +T123+T213
     %     +2.*T123*(PD423/QQ*(QP1/QMP1+1.)+PD123/QQ)
     %     +2.*T213*(PD413/QQ*(QP2/QMP2+1.)+PD213/QQ)
     %     +T341*(PD241/QMP3+1.-2.*PD241/QQ*(QP3/QMP3+1.)
     %           -2.*PD341/QQ)
     %     +T342*(PD142/QMP3+1.-2.*PD142/QQ*(QP3/QMP3+1.)
     %           -2.*PD342/QQ)
 
      BRACK4B=-3./4.*(S2314*(2.*(QP2-QP3)/QQ+1.)
     %             +S1324*(2.*(QP1-QP3)/QQ+1.)
     %             +S1423*(2.*(QP1-QP4)/QQ+1.)
     %             +S2413*(2.*(QP2-QP4)/QQ+1.)
     %             +4.*S34*(QP4-QP3)/QQ)
 
      BRACK4=BRACK4A+BRACK4B
 
      DO 208 K=1,4
 
      HCOMP1(K)=(PIM3(K)-PIM4(K))*BRACK1
      HCOMP2(K)=PIM1(K)*BRACK2
      HCOMP3(K)=PIM2(K)*BRACK3
      HCOMP4(K)=(PIM1(K)+PIM2(K)+PIM3(K)+PIM4(K))*BRACK4
 
 208  CONTINUE
 
      DO 209 I=1,4
 
      HADCUR(I)=HCOMP1(I)-HCOMP2(I)-HCOMP3(I)+HCOMP4(I)
      HADCUR(I)=-COEF1*FRHO4(SQRT(QQ),AMPI,AMPI)*HADCUR(I)
 
 209  CONTINUE
 
 
C --- END OF THE NON OMEGA CURRENT (3 POSSIBILITIES)
 201   CONTINUE
C
C
C --- THERE ARE TWO POSSIBILITIES FOR OMEGA CURRENT
C --- PA PB ARE CORRESPONDING FIRST AND SECOND PI-S
       DO 301 KK=1,2
        DO 302 I=1,4
         PA(I)=PP(KK,I)
         PB(I)=PP(3-KK,I)
 302    CONTINUE
C --- LORENTZ INVARIANTS
         QQA=0.0
         SS23=0.0
         SS24=0.0
         SS34=0.0
         QP1P2=0.0
         QP1P3=0.0
         QP1P4=0.0
         P1P2 =0.0
         P1P3 =0.0
         P1P4 =0.0
        DO 303 K=1,4
                     SIGN=-1.0
         IF (K.EQ.4) SIGN= 1.0
         QQA=QQA+SIGN*(PAA(K)-PA(K))**2
         SS23=SS23+SIGN*(PB(K)  +PIM3(K))**2
         SS24=SS24+SIGN*(PB(K)  +PIM4(K))**2
         SS34=SS34+SIGN*(PIM3(K)+PIM4(K))**2
         QP1P2=QP1P2+SIGN*(PAA(K)-PA(K))*PB(K)
         QP1P3=QP1P3+SIGN*(PAA(K)-PA(K))*PIM3(K)
         QP1P4=QP1P4+SIGN*(PAA(K)-PA(K))*PIM4(K)
         P1P2=P1P2+SIGN*PA(K)*PB(K)
         P1P3=P1P3+SIGN*PA(K)*PIM3(K)
         P1P4=P1P4+SIGN*PA(K)*PIM4(K)
 303    CONTINUE
C
        FORM2=COEF2*(BWIGN(QQ,AMRO,GAMRO)+ELPHA*BWIGN(QQ,AMROP,GAMROP))
C        FORM3=BWIGN(QQA,AMOM,GAMOM)*(BWIGN(SS23,AMRO,GAMRO)+
C     $        BWIGN(SS24,AMRO,GAMRO)+BWIGN(SS34,AMRO,GAMRO))
        FORM3=BWIGN(QQA,AMOM,GAMOM)
C
        DO 304 K=1,4
          HADCUR(K)=HADCUR(K)+FORM2*FORM3*(
     $              PB  (K)*(QP1P3*P1P4-QP1P4*P1P3)
     $             +PIM3(K)*(QP1P4*P1P2-QP1P2*P1P4)
     $             +PIM4(K)*(QP1P2*P1P3-QP1P3*P1P2) )
 304    CONTINUE
 301   CONTINUE
C
      ELSE
C ===================================================================
C PI0 PI0 P0 PI- CASE                                            ====
C ===================================================================
       QQ=PAA(4)**2-PAA(3)**2-PAA(2)**2-PAA(1)**2
 
 
C FIRST DEFINITION OF SCALAR PRODUCTS OF MOMENTUM VECTORS
 
C DEFINE (Q-PI)**2 AS QPI:
 
      QMP1=(PIM2(4)+PIM3(4)+PIM4(4))**2-(PIM2(3)+PIM3(3)+PIM4(3))**2
     %   -(PIM2(2)+PIM3(2)+PIM4(2))**2-(PIM2(1)+PIM3(1)+PIM4(1))**2
 
      QMP2=(PIM1(4)+PIM3(4)+PIM4(4))**2-(PIM1(3)+PIM3(3)+PIM4(3))**2
     %   -(PIM1(2)+PIM3(2)+PIM4(2))**2-(PIM1(1)+PIM3(1)+PIM4(1))**2
 
      QMP3=(PIM1(4)+PIM2(4)+PIM4(4))**2-(PIM1(3)+PIM2(3)+PIM4(3))**2
     %   -(PIM1(2)+PIM2(2)+PIM4(2))**2-(PIM1(1)+PIM2(1)+PIM4(1))**2
 
      QMP4=(PIM1(4)+PIM2(4)+PIM3(4))**2-(PIM1(3)+PIM2(3)+PIM3(3))**2
     %   -(PIM1(2)+PIM2(2)+PIM3(2))**2-(PIM1(1)+PIM2(1)+PIM3(1))**2
 
 
C DEFINE (PI+PK)**2 AS PSIK:
 
      PS14=(PIM1(4)+PIM4(4))**2-(PIM1(3)+PIM4(3))**2
     %    -(PIM1(2)+PIM4(2))**2-(PIM1(1)+PIM4(1))**2
 
      PS21=(PIM2(4)+PIM1(4))**2-(PIM2(3)+PIM1(3))**2
     %    -(PIM2(2)+PIM1(2))**2-(PIM2(1)+PIM1(1))**2
 
      PS23=(PIM2(4)+PIM3(4))**2-(PIM2(3)+PIM3(3))**2
     %    -(PIM2(2)+PIM3(2))**2-(PIM2(1)+PIM3(1))**2
 
      PS24=(PIM2(4)+PIM4(4))**2-(PIM2(3)+PIM4(3))**2
     %    -(PIM2(2)+PIM4(2))**2-(PIM2(1)+PIM4(1))**2
 
      PS31=(PIM3(4)+PIM1(4))**2-(PIM3(3)+PIM1(3))**2
     %    -(PIM3(2)+PIM1(2))**2-(PIM3(1)+PIM1(1))**2
 
      PS34=(PIM3(4)+PIM4(4))**2-(PIM3(3)+PIM4(3))**2
     %    -(PIM3(2)+PIM4(2))**2-(PIM3(1)+PIM4(1))**2
 
 
 
      PD324=PIM3(4)*(PIM2(4)-PIM4(4))-PIM3(3)*(PIM2(3)-PIM4(3))
     %     -PIM3(2)*(PIM2(2)-PIM4(2))-PIM3(1)*(PIM2(1)-PIM4(1))
 
      PD314=PIM3(4)*(PIM1(4)-PIM4(4))-PIM3(3)*(PIM1(3)-PIM4(3))
     %     -PIM3(2)*(PIM1(2)-PIM4(2))-PIM3(1)*(PIM1(1)-PIM4(1))
 
      PD234=PIM2(4)*(PIM3(4)-PIM4(4))-PIM2(3)*(PIM3(3)-PIM4(3))
     %     -PIM2(2)*(PIM3(2)-PIM4(2))-PIM2(1)*(PIM3(1)-PIM4(1))
 
      PD214=PIM2(4)*(PIM1(4)-PIM4(4))-PIM2(3)*(PIM1(3)-PIM4(3))
     %     -PIM2(2)*(PIM1(2)-PIM4(2))-PIM2(1)*(PIM1(1)-PIM4(1))
 
      PD134=PIM1(4)*(PIM3(4)-PIM4(4))-PIM1(3)*(PIM3(3)-PIM4(3))
     %     -PIM1(2)*(PIM3(2)-PIM4(2))-PIM1(1)*(PIM3(1)-PIM4(1))
 
      PD124=PIM1(4)*(PIM2(4)-PIM4(4))-PIM1(3)*(PIM2(3)-PIM4(3))
     %     -PIM1(2)*(PIM2(2)-PIM4(2))-PIM1(1)*(PIM2(1)-PIM4(1))
 
C DEFINE Q*PI = QPI:
 
      QP1=PIM1(4)*(PIM1(4)+PIM2(4)+PIM3(4)+PIM4(4))
     %   -PIM1(3)*(PIM1(3)+PIM2(3)+PIM3(3)+PIM4(3))
     %   -PIM1(2)*(PIM1(2)+PIM2(2)+PIM3(2)+PIM4(2))
     %   -PIM1(1)*(PIM1(1)+PIM2(1)+PIM3(1)+PIM4(1))
 
      QP2=PIM2(4)*(PIM1(4)+PIM2(4)+PIM3(4)+PIM4(4))
     %   -PIM2(3)*(PIM1(3)+PIM2(3)+PIM3(3)+PIM4(3))
     %   -PIM2(2)*(PIM1(2)+PIM2(2)+PIM3(2)+PIM4(2))
     %   -PIM2(1)*(PIM1(1)+PIM2(1)+PIM3(1)+PIM4(1))
 
      QP3=PIM3(4)*(PIM1(4)+PIM2(4)+PIM3(4)+PIM4(4))
     %   -PIM3(3)*(PIM1(3)+PIM2(3)+PIM3(3)+PIM4(3))
     %   -PIM3(2)*(PIM1(2)+PIM2(2)+PIM3(2)+PIM4(2))
     %   -PIM3(1)*(PIM1(1)+PIM2(1)+PIM3(1)+PIM4(1))
 
      QP4=PIM4(4)*(PIM1(4)+PIM2(4)+PIM3(4)+PIM4(4))
     %   -PIM4(3)*(PIM1(3)+PIM2(3)+PIM3(3)+PIM4(3))
     %   -PIM4(2)*(PIM1(2)+PIM2(2)+PIM3(2)+PIM4(2))
     %   -PIM4(1)*(PIM1(1)+PIM2(1)+PIM3(1)+PIM4(1))
 
 
C DEFINE T(PI;PJ,PK)= TIJK:
 
      T324=BWIGA1(QMP3)*FPIKM(SQRT(PS24),AMPI,AMPI)*GAMMA1
      T314=BWIGA1(QMP3)*FPIKM(SQRT(PS14),AMPI,AMPI)*GAMMA1
      T234=BWIGA1(QMP2)*FPIKM(SQRT(PS34),AMPI,AMPI)*GAMMA1
      T214=BWIGA1(QMP2)*FPIKM(SQRT(PS14),AMPI,AMPI)*GAMMA1
      T134=BWIGA1(QMP1)*FPIKM(SQRT(PS34),AMPI,AMPI)*GAMMA1
      T124=BWIGA1(QMP1)*FPIKM(SQRT(PS24),AMPI,AMPI)*GAMMA1
 
C DEFINE S(I,J;K,L)= SIJKL:
 
      S1423=FRHO4(SQRT(PS14),AMPI,AMPI)*BWIGEPS(PS23)*GAMMA2
      S2431=FRHO4(SQRT(PS24),AMPI,AMPI)*BWIGEPS(PS31)*GAMMA2
      S3421=FRHO4(SQRT(PS34),AMPI,AMPI)*BWIGEPS(PS21)*GAMMA2
 
 
C DEFINITION OF AMPLITUDE, FIRST THE [] BRACKETS:
 
      BRACK1=T234+T324+2.*T314+T134+2.*T214+T124
     %    +T134*PD234/QMP1+T124*PD324/QMP1
     %    -3./2.*(S3421+S2431+2.*S1423)
 
 
      BRACK2=T234*(1.+2.*PD134/QMP2)+3.*T324+3.*T124
     %    +T134*(1.-PD234/QMP1)+2.*T214*PD314/QMP2
     %    -T124*PD324/QMP1
     %    -3./2.*(S3421+3.*S2431)
 
      BRACK3=T324*(1.+2.*PD124/QMP3)+3.*T234+3.*T134
     %    +T124*(1.-PD324/QMP1)+2.*T314*PD214/QMP3
     %    -T134*PD234/QMP1
     %    -3./2.*(3.*S3421+S2431)
 
      BRACK4A=2.*T234*(1./2.+PD134/QQ*(QP2/QMP2+1.)+PD234/QQ)
     %     +2.*T324*(1./2.+PD124/QQ*(QP3/QMP3+1.)+PD324/QQ)
     %     +2.*T134*(1./2.+PD234/QQ*(QP1/QMP1+1.)
     %              -1./2.*PD234/QMP1+PD134/QQ)
     %     +2.*T124*(1./2.+PD324/QQ*(QP1/QMP1+1.)
     %              -1./2.*PD324/QMP1+PD124/QQ)
     %     +2.*T214*(PD314/QQ*(QP2/QMP2+1.)+PD214/QQ)
     %     +2.*T314*(PD214/QQ*(QP3/QMP3+1.)+PD314/QQ)
 
      BRACK4B=-3./2.*(S3421*(2.*(QP3-QP4)/QQ+1.)
     %             +S2431*(2.*(QP2-QP4)/QQ+1.)
     %             +S1423*2.*(QP1-QP4)/QQ)
 
 
      BRACK4=BRACK4A+BRACK4B
 
      DO 308 K=1,4
 
      HCOMP1(K)=(PIM1(K)-PIM4(K))*BRACK1
      HCOMP2(K)=PIM2(K)*BRACK2
      HCOMP3(K)=PIM3(K)*BRACK3
      HCOMP4(K)=(PIM1(K)+PIM2(K)+PIM3(K)+PIM4(K))*BRACK4
 
 308  CONTINUE
 
      DO 309 I=1,4
 
      HADCUR(I)=HCOMP1(I)+HCOMP2(I)+HCOMP3(I)-HCOMP4(I)
      HADCUR(I)=COEF1*FRHO4(SQRT(QQ),AMPI,AMPI)*HADCUR(I)
 
 309  CONTINUE
 
 101   CONTINUE
      ENDIF
C M. Finkemeier et al. END
      END
