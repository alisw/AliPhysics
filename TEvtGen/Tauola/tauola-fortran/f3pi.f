

*AJW 1 version of a1 form factor
      COMPLEX FUNCTION F3PI(IFORM,QQ,SA,SB)
C.......................................................................
C.
C. F3PI - 1 version of a1 form factor, used in TAUOLA
C.
C. Inputs    : None
C.           :
C. Outputs   : None
C.
C. COMMON    : None
C.
C. Calls     : 
C. Called    : by FORM1-FORM3 in $C_CVSSRC/korb/koralb/formf.F
C. Author    : Alan Weinstein 2/98
C.
C. Detailed description
C.   First determine whether we are doing pi-2pi0 or 3pi.
C.   Then implement full form-factor from fit:
C.   [(rho-pi S-wave) + (rho-prim-pi S-wave) +
C.    (rho-pi D-wave) + (rho-prim-pi D-wave) + 
C.    (f2 pi D-wave) + (sigmapi S-wave) + (f0pi S-wave)]
C.   based on fit to pi-2pi0 by M. Schmidler, CBX 97-64-Update (4/22/98)
C.   All the parameters in this routine are hard-coded!!
C.
C.......................................................................



* -------------------- Argument declarations ---------------
 
      INTEGER IFORM
      REAL QQ,SA,SB
* -------------------- EXTERNAL declarations ---------------
*
      REAL PKORB
      COMPLEX BWIGML
* -------------------- SEQUENCE declarations ---------------
*
* -------------------- Local    declarations ---------------
*
      CHARACTER*(*) CRNAME
      PARAMETER(    CRNAME = 'F3PI' )
*
      INTEGER IFIRST,IDK
      REAL MRO,GRO,MRP,GRP,MF2,GF2,MF0,GF0,MSG,GSG
      REAL M1,M2,M3,M1SQ,M2SQ,M3SQ,MPIZ,MPIC
      REAL S1,S2,S3,R,PI
      REAL F134,F150,F15A,F15B,F167
      REAL F34A,F34B,F35,F35A,F35B,F36A,F36B
      COMPLEX BT1,BT2,BT3,BT4,BT5,BT6,BT7
      COMPLEX FRO1,FRO2,FRP1,FRP2
      COMPLEX FF21,FF22,FF23,FSG1,FSG2,FSG3,FF01,FF02,FF03
      COMPLEX FA1A1P,FORMA1

* -------------------- SAVE     declarations ---------------
*
* -------------------- DATA  initializations ---------------
*
      DATA IFIRST/0/
* ----------------- Executable code starts here ------------
*
C. Hard-code the fit parameters:
      IF (IFIRST.EQ.0) THEN
        IFIRST = 1
C rho, rhoprime, f2(1275), f0(1186), sigma(made up!)
        MRO = 0.7743
        GRO = 0.1491
        MRP = 1.370 
        GRP = 0.386 
        MF2 = 1.275
        GF2 = 0.185
        MF0 = 1.186
        GF0 = 0.350
        MSG = 0.860
        GSG = 0.880
        MPIZ = PKORB(1,7)
        MPIC = PKORB(1,8)

C Fit coefficients for each of the contributions:
        PI = 3.14159
        BT1 = CMPLX(1.,0.)
        BT2 = CMPLX(0.12,0.)*CEXP(CMPLX(0., 0.99*PI))
        BT3 = CMPLX(0.37,0.)*CEXP(CMPLX(0.,-0.15*PI))
        BT4 = CMPLX(0.87,0.)*CEXP(CMPLX(0., 0.53*PI))
        BT5 = CMPLX(0.71,0.)*CEXP(CMPLX(0., 0.56*PI))
        BT6 = CMPLX(2.10,0.)*CEXP(CMPLX(0., 0.23*PI))
        BT7 = CMPLX(0.77,0.)*CEXP(CMPLX(0.,-0.54*PI))

        PRINT *,' In F3pi: add (rho-pi S-wave) + (rhop-pi S-wave) +'
        PRINT *,'              (rho-pi D-wave) + (rhop-pi D-wave) +'
        PRINT *,'   (f2 pi D-wave) + (sigmapi S-wave) + (f0pi S-wave)'
      END IF

C Initialize to 0:
      F3PI = CMPLX(0.,0.)

C.   First determine whether we are doing pi-2pi0 or 3pi.
C     PKORB is set up to remember what flavor of 3pi it gave to KORALB,
C     since KORALB doesnt bother to remember!!
      R = PKORB(4,11)
      IF (R.EQ.0.) THEN
C it is 2pi0pi-
        IDK = 1
        M1 = MPIZ
        M2 = MPIZ
        M3 = MPIC
      ELSE
C it is 3pi
        IDK = 2
        M1 = MPIC
        M2 = MPIC
        M3 = MPIC
      END IF
      M1SQ = M1*M1
      M2SQ = M2*M2
      M3SQ = M3*M3

C.   Then implement full form-factor from fit:
C.   [(rho-pi S-wave) + (rho-prim-pi S-wave) +
C.    (rho-pi D-wave) + (rho-prim-pi D-wave) + 
C.    (f2 pi D-wave) + (sigmapi S-wave) + (f0pi S-wave)]
C.   based on fit to pi-2pi0 by M. Schmidler, CBX 97-64-Update (4/22/98)

C Note that for FORM1, the arguments are S1, S2;
C           for FORM2, the arguments are S2, S1;
C           for FORM3, the arguments are S3, S1.
C Here, we implement FORM1 and FORM2 at the same time,
C so the above switch is just what we need!

      IF (IFORM.EQ.1.OR.IFORM.EQ.2) THEN
        S1 = SA
        S2 = SB
        S3 = QQ-SA-SB+M1SQ+M2SQ+M3SQ
        IF (S3.LE.0..OR.S2.LE.0.) RETURN

        IF (IDK.EQ.1) THEN
C it is 2pi0pi-
C Lorentz invariants for all the contributions:
          F134 = -(1./3.)*((S3-M3SQ)-(S1-M1SQ))
          F150 =  (1./18.)*(QQ-M3SQ+S3)*(2.*M1SQ+2.*M2SQ-S3)/S3
          F167 =  (2./3.)

C Breit Wigners for all the contributions:
          FRO1 = BWIGML(S1,MRO,GRO,M2,M3,1)
          FRP1 = BWIGML(S1,MRP,GRP,M2,M3,1)
          FRO2 = BWIGML(S2,MRO,GRO,M3,M1,1)
          FRP2 = BWIGML(S2,MRP,GRP,M3,M1,1)
          FF23 = BWIGML(S3,MF2,GF2,M1,M2,2)
          FSG3 = BWIGML(S3,MSG,GSG,M1,M2,0)
          FF03 = BWIGML(S3,MF0,GF0,M1,M2,0)

          F3PI = BT1*FRO1+BT2*FRP1+
     1       BT3*CMPLX(F134,0.)*FRO2+BT4*CMPLX(F134,0.)*FRP2+
     1       BT5*CMPLX(F150,0.)*FF23+
     1       BT6*CMPLX(F167,0.)*FSG3+BT7*CMPLX(F167,0.)*FF03

C         F3PI = FPIKM(SQRT(S1),M2,M3)
        ELSEIF (IDK.EQ.2) THEN
C it is 3pi
C Lorentz invariants for all the contributions:
          F134 = -(1./3.)*((S3-M3SQ)-(S1-M1SQ))
          F15A = -(1./2.)*((S2-M2SQ)-(S3-M3SQ))
          F15B = -(1./18.)*(QQ-M2SQ+S2)*(2.*M1SQ+2.*M3SQ-S2)/S2
          F167 = -(2./3.)

C Breit Wigners for all the contributions:
          FRO1 = BWIGML(S1,MRO,GRO,M2,M3,1)
          FRP1 = BWIGML(S1,MRP,GRP,M2,M3,1)
          FRO2 = BWIGML(S2,MRO,GRO,M3,M1,1)
          FRP2 = BWIGML(S2,MRP,GRP,M3,M1,1)
          FF21 = BWIGML(S1,MF2,GF2,M2,M3,2)
          FF22 = BWIGML(S2,MF2,GF2,M3,M1,2)
          FSG2 = BWIGML(S2,MSG,GSG,M3,M1,0)
          FF02 = BWIGML(S2,MF0,GF0,M3,M1,0)

          F3PI = BT1*FRO1+BT2*FRP1+
     1       BT3*CMPLX(F134,0.)*FRO2+BT4*CMPLX(F134,0.)*FRP2
     1      -BT5*CMPLX(F15A,0.)*FF21-BT5*CMPLX(F15B,0.)*FF22
     1      -BT6*CMPLX(F167,0.)*FSG2-BT7*CMPLX(F167,0.)*FF02

C         F3PI = FPIKM(SQRT(S1),M2,M3)
        END IF

      ELSE IF (IFORM.EQ.3) THEN
        S3 = SA
        S1 = SB
        S2 = QQ-SA-SB+M1SQ+M2SQ+M3SQ
        IF (S1.LE.0..OR.S2.LE.0.) RETURN

        IF (IDK.EQ.1) THEN
C it is 2pi0pi-
C Lorentz invariants for all the contributions:
          F34A = (1./3.)*((S2-M2SQ)-(S3-M3SQ))
          F34B = (1./3.)*((S3-M3SQ)-(S1-M1SQ))
          F35  =-(1./2.)*((S1-M1SQ)-(S2-M2SQ))

C Breit Wigners for all the contributions:
          FRO1 = BWIGML(S1,MRO,GRO,M2,M3,1)
          FRP1 = BWIGML(S1,MRP,GRP,M2,M3,1)
          FRO2 = BWIGML(S2,MRO,GRO,M3,M1,1)
          FRP2 = BWIGML(S2,MRP,GRP,M3,M1,1)
          FF23 = BWIGML(S3,MF2,GF2,M1,M2,2)

          F3PI = 
     1       BT3*(CMPLX(F34A,0.)*FRO1+CMPLX(F34B,0.)*FRO2)+
     1       BT4*(CMPLX(F34A,0.)*FRP1+CMPLX(F34B,0.)*FRP2)+
     1       BT5*CMPLX(F35,0.)*FF23
     
C         F3PI = CMPLX(0.,0.)
        ELSEIF (IDK.EQ.2) THEN
C it is 3pi
C Lorentz invariants for all the contributions:
          F34A = (1./3.)*((S2-M2SQ)-(S3-M3SQ))
          F34B = (1./3.)*((S3-M3SQ)-(S1-M1SQ))
          F35A = -(1./18.)*(QQ-M1SQ+S1)*(2.*M2SQ+2.*M3SQ-S1)/S1
          F35B =  (1./18.)*(QQ-M2SQ+S2)*(2.*M3SQ+2.*M1SQ-S2)/S2
          F36A = -(2./3.)
          F36B =  (2./3.)

C Breit Wigners for all the contributions:
          FRO1 = BWIGML(S1,MRO,GRO,M2,M3,1)
          FRP1 = BWIGML(S1,MRP,GRP,M2,M3,1)
          FRO2 = BWIGML(S2,MRO,GRO,M3,M1,1)
          FRP2 = BWIGML(S2,MRP,GRP,M3,M1,1)
          FF21 = BWIGML(S1,MF2,GF2,M2,M3,2)
          FF22 = BWIGML(S2,MF2,GF2,M3,M1,2)
          FSG1 = BWIGML(S1,MSG,GSG,M2,M3,0)
          FSG2 = BWIGML(S2,MSG,GSG,M3,M1,0)
          FF01 = BWIGML(S1,MF0,GF0,M2,M3,0)
          FF02 = BWIGML(S2,MF0,GF0,M3,M1,0)

          F3PI = 
     1       BT3*(CMPLX(F34A,0.)*FRO1+CMPLX(F34B,0.)*FRO2)+
     1       BT4*(CMPLX(F34A,0.)*FRP1+CMPLX(F34B,0.)*FRP2)
     1      -BT5*(CMPLX(F35A,0.)*FF21+CMPLX(F35B,0.)*FF22)
     1      -BT6*(CMPLX(F36A,0.)*FSG1+CMPLX(F36B,0.)*FSG2)
     1      -BT7*(CMPLX(F36A,0.)*FF01+CMPLX(F36B,0.)*FF02)
     
C         F3PI = CMPLX(0.,0.)
        END IF
      END IF

C Add overall a1/a1prime:
      FORMA1 = FA1A1P(QQ)
      F3PI = F3PI*FORMA1

      RETURN
      END
C **********************************************************
      COMPLEX FUNCTION BWIGML(S,M,G,M1,M2,L)
C **********************************************************
C     L-WAVE BREIT-WIGNER
C **********************************************************
      REAL S,M,G,M1,M2
      INTEGER L,IPOW
      REAL MSQ,W,WGS,MP,MM,QS,QM

      MP = (M1+M2)**2
      MM = (M1-M2)**2
      MSQ = M*M
      W = SQRT(S)
      WGS = 0.0
      IF (W.GT.(M1+M2)) THEN
        QS=SQRT(ABS((S   -MP)*(S   -MM)))/W
        QM=SQRT(ABS((MSQ -MP)*(MSQ -MM)))/M
        IPOW = 2*L+1
        WGS=G*(MSQ/W)*(QS/QM)**IPOW
      ENDIF

      BWIGML=CMPLX(MSQ,0.)/CMPLX(MSQ-S,-WGS)

      RETURN
      END
C=======================================================================
      COMPLEX FUNCTION FA1A1P(XMSQ)
C     ==================================================================
C     complex form-factor for a1+a1prime.                       AJW 1/98
C     ==================================================================

      REAL XMSQ
      REAL PKORB,WGA1
      REAL XM1,XG1,XM2,XG2,XM1SQ,XM2SQ,GG1,GG2,GF,FG1,FG2
      COMPLEX BET,F1,F2
      INTEGER IFIRST/0/

      IF (IFIRST.EQ.0) THEN
        IFIRST = 1

C The user may choose masses and widths that differ from nominal:
        XM1 = PKORB(1,10)
        XG1 = PKORB(2,10)
        XM2 = PKORB(1,17)
        XG2 = PKORB(2,17)
        BET = CMPLX(PKORB(3,17),0.)
C scale factors relative to nominal:
        GG1 = XM1*XG1/(1.3281*0.806)
        GG2 = XM2*XG2/(1.3281*0.806)

        XM1SQ = XM1*XM1
        XM2SQ = XM2*XM2
      END IF

      GF = WGA1(XMSQ)
      FG1 = GG1*GF
      FG2 = GG2*GF
      F1 = CMPLX(-XM1SQ,0.0)/CMPLX(XMSQ-XM1SQ,FG1)
      F2 = CMPLX(-XM2SQ,0.0)/CMPLX(XMSQ-XM2SQ,FG2)
      FA1A1P = F1+BET*F2

      RETURN
      END
C=======================================================================
      FUNCTION WGA1(QQ)

C mass-dependent M*Gamma of a1 through its decays to 
C.   [(rho-pi S-wave) + (rho-pi D-wave) + 
C.    (f2 pi D-wave) + (f0pi S-wave)]
C.  AND simple K*K S-wave

      REAL QQ,WGA1
      DOUBLE PRECISION MKST,MK,MK1SQ,MK2SQ,C3PI,CKST
      DOUBLE PRECISION S,WGA1C,WGA1N,WG3PIC,WG3PIN,GKST
      INTEGER IFIRST
C-----------------------------------------------------------------------
C
      IF (IFIRST.NE.987) THEN
        IFIRST = 987
C
C Contribution to M*Gamma(m(3pi)^2) from S-wave K*K:
        MKST = 0.894D0
        MK = 0.496D0
        MK1SQ = (MKST+MK)**2
        MK2SQ = (MKST-MK)**2
C coupling constants squared:
        C3PI = 0.2384D0**2
        CKST = 4.7621D0**2*C3PI
      END IF

C-----------------------------------------------------------------------
C Parameterization of numerical integral of total width of a1 to 3pi.
C From M. Schmidtler, CBX-97-64-Update.
      S = DBLE(QQ)
      WG3PIC = WGA1C(S)
      WG3PIN = WGA1N(S)

C Contribution to M*Gamma(m(3pi)^2) from S-wave K*K, if above threshold
      GKST = 0.D0
      IF (S.GT.MK1SQ) GKST = SQRT((S-MK1SQ)*(S-MK2SQ))/(2.*S)

      WGA1 = SNGL(C3PI*(WG3PIC+WG3PIN)+CKST*GKST)

      RETURN 
      END
C=======================================================================
      DOUBLE PRECISION FUNCTION WGA1C(S)
C
C parameterization of m*Gamma(m^2) for pi-2pi0 system
C
      DOUBLE PRECISION S,STH,Q0,Q1,Q2,P0,P1,P2,P3,P4,G1_IM
C
      PARAMETER(Q0 =   5.80900D0,Q1 =  -3.00980D0,Q2 =   4.57920D0,
     1          P0 = -13.91400D0,P1 =  27.67900D0,P2 = -13.39300D0,
     2          P3 =   3.19240D0,P4 =  -0.10487D0)
C
      PARAMETER (STH   = 0.1753D0)
C---------------------------------------------------------------------

      IF(S.LT.STH) THEN
       G1_IM = 0.D0
      ELSEIF((S.GT.STH).AND.(S.LT.0.823D0)) THEN
       G1_IM = Q0*(S-STH)**3*(1. + Q1*(S-STH) + Q2*(S-STH)**2)
      ELSE
       G1_IM = P0 + P1*S + P2*S**2+ P3*S**3 + P4*S**4
      ENDIF

      WGA1C = G1_IM      
      RETURN
      END
C=======================================================================
      DOUBLE PRECISION FUNCTION WGA1N(S)
C
C parameterization of m*Gamma(m^2) for pi-pi+pi- system
C
      DOUBLE PRECISION S,STH,Q0,Q1,Q2,P0,P1,P2,P3,P4,G1_IM
C
      PARAMETER(Q0 =   6.28450D0,Q1 =  -2.95950D0,Q2 =   4.33550D0,
     1          P0 = -15.41100D0,P1 =  32.08800D0,P2 = -17.66600D0,
     2          P3 =   4.93550D0,P4 =  -0.37498D0)
C
      PARAMETER (STH   = 0.1676D0)
C---------------------------------------------------------------------

      IF(S.LT.STH) THEN
       G1_IM = 0.D0
      ELSEIF((S.GT.STH).AND.(S.LT.0.823D0)) THEN
       G1_IM = Q0*(S-STH)**3*(1. + Q1*(S-STH) + Q2*(S-STH)**2)
      ELSE
       G1_IM = P0 + P1*S + P2*S**2+ P3*S**3 + P4*S**4
      ENDIF

      WGA1N = G1_IM      
      RETURN
      END
