*CMZ :          17/07/98  15.49.05  by  Federico Carminati
*-- Author :
      SUBROUTINE SHPROD
c	=================

c	Phase Space Particle Generation

*KEEP,SHPHYP.
      COMMON /SHPHYP/ JWEI,NDNDY,YLIM,PTLIM,JWEAK,JPI0,JETA,JPIC,JPRO,
     +                  JKAC,JKA0,JRHO,JOME,JPHI,JPSI,JDRY
*KEEP,SHGENE.
      COMMON /SHGENE/ IEVT,NPI0,NETA,NPIC,NPRO,NKAC,NKA0,NRHO,NOME,
     +                  NPHI,NPSI,NDRY
*KEEP,SHPRAT.
      COMMON /SHPRAT/ PI0R,ETAR,RHOR,OMER,PHIR,PSIR,DRYR
*KEEP,SHNORM.
      COMMON /SHNORM/ PINOR,PIRAT,ETANOR,ETARAT,RHONOR,OMENOR,PHINOR,
     +                  PSINOR,DRYNOR
*KEEP,LUJETS.
      COMMON /LUJETS/ N,K(200000,5),P(200000,5),V(200000,5)
      SAVE /LUJETS/
*KEEP,LUDAT1.
      COMMON /LUDAT1/ MSTU(200),PARU(200),MSTJ(200),PARJ(200)
      SAVE /LUDAT1/
*KEND.

      DIMENSION PDY(4),PDA(4),PDB(4)

      N = 0

c	pi0 generation

      DO 10 IPI0=1,NPI0
        CALL SHRAPI(Y)
        CALL SHPTGE(111,PT,PHI,W)
        CALL SHTOLU(111,Y,PT,PHI,W)
10	CONTINUE

c	eta generation

      DO 20 IETA=1,NETA
        CALL SHRAPI(Y)
        CALL SHPTGE(221,PT,PHI,W)
        CALL SHTOLU(221,Y,PT,PHI,W)
20	CONTINUE

c	pion generation

      DO 30 IPIC=1,NPIC
        CALL SHRAPI(Y)
        CALL SHPTGE(211,PT,PHI,W)
        IF (RLU(0).LE..5) THEN
          CALL SHTOLU(-211,Y,PT,PHI,W)
        ELSE
          CALL SHTOLU(211,Y,PT,PHI,W)
        ENDIF
30	CONTINUE

c	proton generation

      DO 40 IPRO=1,NPRO
        CALL SHRAPI(Y)
        CALL SHPTGE(2212,PT,PHI,W)
        IF (RLU(0).LE..5) THEN
          CALL SHTOLU(-2212,Y,PT,PHI,W)
        ELSE
          CALL SHTOLU(2212,Y,PT,PHI,W)
        ENDIF
40	CONTINUE

c	kaon generation

      DO 50 IKAC=1,NKAC
        CALL SHRAPI(Y)
        CALL SHPTGE(321,PT,PHI,W)
        IF (RLU(0).LE..5) THEN
          CALL SHTOLU(-321,Y,PT,PHI,W)
        ELSE
          CALL SHTOLU(321,Y,PT,PHI,W)
        ENDIF
50	CONTINUE

c	K0 generation

      DO 60 IKA0=1,NKA0
        CALL SHRAPI(Y)
        CALL SHPTGE(311,PT,PHI,W)
        IF (RLU(0).LE..5) THEN
          CALL SHTOLU(-311,Y,PT,PHI,W)
        ELSE
          CALL SHTOLU(311,Y,PT,PHI,W)
        ENDIF
60	CONTINUE

c	rho generation

      DO 70 IRHO=1,NRHO
        CALL SHRAPI(Y)
        CALL SHPTGE(113,PT,PHI,W)
        CALL SHTOLU(113,Y,PT,PHI,W)
70	CONTINUE

c	omega generation

      DO 80 IOME=1,NOME
        CALL SHRAPI(Y)
        CALL SHPTGE(223,PT,PHI,W)
        CALL SHTOLU(223,Y,PT,PHI,W)
80	CONTINUE

c	phi generation

      DO 90 IPHI=1,NPHI
        CALL SHRAPI(Y)
        CALL SHPTGE(333,PT,PHI,W)
        CALL SHTOLU(333,Y,PT,PHI,W)
90	CONTINUE

c	J/psi generation

      DO 100 IPSI=1,NPSI
        CALL SHRAPI(Y)
        CALL SHPTGE(443,PT,PHI,W)
        CALL SHTOLU(443,Y,PT,PHI,W)
100	CONTINUE

c	Drell-Yan generation

      DO 110 IDRY=1,NDRY

        CALL SHRAPI(Y)
        AMDY2 = 1/(1-RLU(0))
        AMDY  = SQRT(AMDY2)

        PHI = 3.14159*2*(RLU(0)-.5)

        IF (JWEI.EQ.1) THEN
          PT = PTLIM*RLU(0)
          W = PTLIM*SHMTSC(AMDY,PT)*DRYR/DRYNOR/FLOAT(NDRY)
        ELSE
            WRITE(MSTU(11),*) 'ERROR:'
            WRITE(MSTU(11),*) 'Drell-Yan NOT generated with JWEI=1!'
          WRITE(MSTU(11),*)'EXECUTION STOPPED!'
          STOP
        ENDIF

        TM = SQRT(AMDY2+PT**2)
        PDY(1) = PT*COS(PHI)
        PDY(2) = PT*SIN(PHI)
        PDY(3) = TM*SINH(Y)
        PDY(4) = TM*COSH(Y)

        AMDA = ULMASS(11)
        AMDB = ULMASS(11)
        CALL SH2BOD(AMDY,AMDA,AMDB,PDY,PDA,PDB)

        EDA = SQRT(VDOT(PDA,PDA,3)+AMDA**2)
        YDA  = .5 * LOG ( (EDA+PDA(3)) / (EDA-PDA(3)) )
        PTDA = SQRT(VDOT(PDA,PDA,2))
        CA = PDA(1)/PTDA
        SA = PDA(2)/PTDA
        PHIA = ATAN2(SA,CA)
        CALL SHTOLU(-11,YDA,PTDA,PHIA,W)

        EDB = SQRT(VDOT(PDB,PDB,3)+AMDB**2)
        YDB  = .5 * LOG ( (EDB+PDB(3)) / (EDB-PDB(3)) )
        PTDB = SQRT(VDOT(PDB,PDB,2))
        CB = PDB(1)/PTDB
        SB = PDB(2)/PTDB
        PHIB = ATAN2(SB,CB)
        CALL SHTOLU(11,YDB,PTDB,PHIB,W)

110	CONTINUE


      RETURN
      END
