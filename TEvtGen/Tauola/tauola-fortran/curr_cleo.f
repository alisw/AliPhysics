

*AJW 1 version of CURR from KORALB.
      SUBROUTINE CURR_CLEO(MNUM,PIM1,PIM2,PIM3,PIM4,HADCUR)
C     ==================================================================
C AJW, 11/97 - based on original CURR from TAUOLA:
C     hadronic current for 4 pi final state
C     R. Fisher, J. Wess and F. Wagner Z. Phys C3 (1980) 313
C     R. Decker Z. Phys C36 (1987) 487.
C     M. Gell-Mann, D. Sharp, W. Wagner Phys. Rev. Lett 8 (1962) 261.
C BUT, rewritten to be more general and less "theoretical",
C  using parameters tuned by Vasia and DSC.
C     ==================================================================
 
      COMMON / PARMAS / AMTAU,AMNUTA,AMEL,AMNUE,AMMU,AMNUMU
     *                 ,AMPIZ,AMPI,AMRO,GAMRO,AMA1,GAMA1
     *                 ,AMK,AMKZ,AMKST,GAMKST
C
      REAL*4            AMTAU,AMNUTA,AMEL,AMNUE,AMMU,AMNUMU
     *                 ,AMPIZ,AMPI,AMRO,GAMRO,AMA1,GAMA1
     *                 ,AMK,AMKZ,AMKST,GAMKST
C
      REAL  PIM1(4),PIM2(4),PIM3(4),PIM4(4)
      COMPLEX HADCUR(4)

      INTEGER K,L,MNUM,K1,K2,IRO,I,J,KK
      REAL PA(4),PB(4),PAA(4)
      REAL AA(4,4),PP(4,4)
      REAL A,XM,XG,G1,G2,G,AMRO2,GAMRO2,AMRO3,GAMRO3,AMOM,GAMOM
      REAL FRO,COEF1,FPI,COEF2,QQ,SK,DENOM,SIG,QQA,SS23,SS24,SS34,QP1P2
      REAL QP1P3,QP1P4,P1P2,P1P3,P1P4,SIGN
      REAL PKORB,AMPA
      COMPLEX ALF0,ALF1,ALF2,ALF3
      COMPLEX LAM0,LAM1,LAM2,LAM3
      COMPLEX BET1,BET2,BET3
      COMPLEX FORM1,FORM2,FORM3,FORM4,FORM2PI
      COMPLEX BWIGM,WIGFOR,FPIKM,FPIKMD
      COMPLEX AMPL(7),AMPR
      COMPLEX BWIGN
C
      BWIGN(A,XM,XG)=1.0/CMPLX(A-XM**2,XM*XG)
C*******************************************************************************
C
C --- masses and constants
      IF (G1.NE.12.924) THEN
      G1=12.924
      G2=1475.98
      FPI=93.3E-3
      G =G1*G2
      FRO=0.266*AMRO**2
      COEF1=2.0*SQRT(3.0)/FPI**2
      COEF2=FRO*G ! overall constant for the omega current
      COEF2= COEF2*0.56  ! factor 0.56 reduces contribution of omega from 68% to 40 %

C masses and widths for for rho-prim and rho-bis:
      AMRO2 = 1.465
      GAMRO2= 0.310
      AMRO3=1.700
      GAMRO3=0.235
C
      AMOM  = PKORB(1,14)
      GAMOM = PKORB(2,14)
      AMRO2 = PKORB(1,21)
      GAMRO2= PKORB(2,21)
      AMRO3 = PKORB(1,22)
      GAMRO3= PKORB(2,22)
C
C Amplitudes for (pi-pi-pi0pi+) -> PS, rho0, rho-, rho+, omega.
      AMPL(1) = CMPLX(PKORB(3,31)*COEF1,0.)
      AMPL(2) = CMPLX(PKORB(3,32)*COEF1,0.)*CEXP(CMPLX(0.,PKORB(3,42)))
      AMPL(3) = CMPLX(PKORB(3,33)*COEF1,0.)*CEXP(CMPLX(0.,PKORB(3,43)))
      AMPL(4) = CMPLX(PKORB(3,34)*COEF1,0.)*CEXP(CMPLX(0.,PKORB(3,44)))
      AMPL(5) = CMPLX(PKORB(3,35)*COEF2,0.)*CEXP(CMPLX(0.,PKORB(3,45)))
C Amplitudes for (pi0pi0pi0pi-) -> PS, rho-.
      AMPL(6) = CMPLX(PKORB(3,36)*COEF1)
      AMPL(7) = CMPLX(PKORB(3,37)*COEF1)
C
C rho' contributions to rho' -> pi-omega:
      ALF0 = CMPLX(PKORB(3,51),0.0)
      ALF1 = CMPLX(PKORB(3,52)*AMRO**2,0.0)
      ALF2 = CMPLX(PKORB(3,53)*AMRO2**2,0.0)
      ALF3 = CMPLX(PKORB(3,54)*AMRO3**2,0.0)
C rho' contribtions to rho' -> rhopipi:
      LAM0 = CMPLX(PKORB(3,55),0.0)
      LAM1 = CMPLX(PKORB(3,56)*AMRO**2,0.0)
      LAM2 = CMPLX(PKORB(3,57)*AMRO2**2,0.0)
      LAM3 = CMPLX(PKORB(3,58)*AMRO3**2,0.0)
C rho contributions to rhopipi, rho -> 2pi:
      BET1 = CMPLX(PKORB(3,59)*AMRO**2,0.0)
      BET2 = CMPLX(PKORB(3,60)*AMRO2**2,0.0)
      BET3 = CMPLX(PKORB(3,61)*AMRO3**2,0.0)
C
      END IF
C**************************************************
C
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

C  Add M(4pi)-dependence to rhopipi channels:
       FORM4= LAM0+LAM1*BWIGN(QQ,AMRO,GAMRO)
     *            +LAM2*BWIGN(QQ,AMRO2,GAMRO2)
     *            +LAM3*BWIGN(QQ,AMRO3,GAMRO3)

C --- loop over five contributions of the rho-pi-pi
       DO 201 K1=1,3
       DO 201 K2=3,4
C
         IF (K2.EQ.K1) THEN
           GOTO 201
         ELSEIF (K2.EQ.3) THEN
C rho-
            AMPR = AMPL(3)
            AMPA = AMPIZ
         ELSEIF (K1.EQ.3) THEN
C rho+
            AMPR = AMPL(4)
            AMPA = AMPIZ
         ELSE
C rho0
            AMPR = AMPL(2)
            AMPA = AMPI
         END IF
C
         SK=(PP(K1,4)+PP(K2,4))**2-(PP(K1,3)+PP(K2,3))**2
     $     -(PP(K1,2)+PP(K2,2))**2-(PP(K1,1)+PP(K2,1))**2

C -- definition of AA matrix
C -- cronecker delta
        DO 202 I=1,4
         DO 203 J=1,4
 203     AA(I,J)=0.0
 202    AA(I,I)=1.0
C ... and the rest ...
        DO 204 L=1,4
         IF (L.NE.K1.AND.L.NE.K2) THEN
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
C
C --- lets add something to HADCURR
C        FORM1= FPIKM(SQRT(SK),AMPI,AMPI) *FPIKM(SQRT(QQ),AMPI,AMPI)
C        FORM1= AMPL(1)+AMPR*FPIKM(SQRT(SK),AMPI,AMPI)

        FORM2PI= BET1*BWIGM(SK,AMRO,GAMRO,AMPA,AMPI)
     1          +BET2*BWIGM(SK,AMRO2,GAMRO2,AMPA,AMPI)
     2          +BET3*BWIGM(SK,AMRO3,GAMRO3,AMPA,AMPI)
        FORM1= AMPL(1)+AMPR*FORM2PI
C
       DO 206 I=1,4
       DO 206 J=1,4
        HADCUR(I)=HADCUR(I)+FORM1*FORM4*AA(I,J)*(PP(K1,J)-PP(K2,J))
 206   CONTINUE
C --- end of the rho-pi-pi current (5 possibilities)
 201   CONTINUE
C
C ===================================================================
C Now modify the coefficient for the omega-pi current:              =
C ===================================================================
       IF (AMPL(5).EQ.CMPLX(0.,0.)) GOTO 311

C Overall rho+rhoprime for the 4pi system:
C       FORM2=AMPL(5)*(BWIGN(QQ,AMRO,GAMRO)+ELPHA*BWIGN(QQ,AMROP,GAMROP))
C Modified M(4pi)-dependence:
       FORM2=AMPL(5)*(ALF0+ALF1*BWIGN(QQ,AMRO,GAMRO)
     *                    +ALF2*BWIGN(QQ,AMRO2,GAMRO2)
     *                    +ALF3*BWIGN(QQ,AMRO3,GAMRO3))
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
C omega -> rho pi for the 3pi system:
C       FORM3=BWIGN(QQA,AMOM,GAMOM)*(BWIGN(SS23,AMRO,GAMRO)+
C     $        BWIGN(SS24,AMRO,GAMRO)+BWIGN(SS34,AMRO,GAMRO))
C No omega -> rho pi; just straight omega:
        FORM3=BWIGN(QQA,AMOM,GAMOM)
C
        DO 304 K=1,4
         HADCUR(K)=HADCUR(K)+FORM2*FORM3*(
     $             PB  (K)*(QP1P3*P1P4-QP1P4*P1P3)
     $            +PIM3(K)*(QP1P4*P1P2-QP1P2*P1P4)
     $            +PIM4(K)*(QP1P2*P1P3-QP1P3*P1P2) )
 304    CONTINUE
 301   CONTINUE
 311   CONTINUE
C
      ELSE
C ===================================================================
C pi0 pi0 p0 pi- case                                            ====
C ===================================================================
       QQ=PAA(4)**2-PAA(3)**2-PAA(2)**2-PAA(1)**2

C --- loop over three contribution of the non-omega current
       DO 101 K=1,3
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
C       FORM1= FPIKM(SQRT(SK),AMPI,AMPI) *FPIKMD(SQRT(QQ),AMPI,AMPI)
CCCCCCCCCCCCC       FORM1=WIGFOR(SK,AMRO,GAMRO)        (tests)
C       FORM1= FPIKM(SQRT(SK),AMPI,AMPI) *FPIKM(SQRT(QQ),AMPI,AMPI)
       FORM1 = AMPL(6)+AMPL(7)*FPIKM(SQRT(SK),AMPI,AMPI)

        DO 106 I=1,4
        DO 106 J=1,4
         HADCUR(I)=HADCUR(I)+FORM1*AA(I,J)*(PP(K,J)-PP(4,J))
 106    CONTINUE
C --- end of the non omega current (3 possibilities)
 101   CONTINUE

      ENDIF
      END
 


