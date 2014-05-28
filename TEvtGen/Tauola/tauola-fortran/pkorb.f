

      REAL FUNCTION PKORB(IF1,IF2)
**********************************************************************
*
* This function returns a real value
* needed in the 1 version of KORALB/TAUOLA
* corresponding to a mass, width, mixing amplitude, or branching fraction
* depending on whether IF1 = 1, 2, 3, 4 respectively.
* The idea is to make minimal mods to the 3-rd party KORALB/TAUOLA code,
* so this function supplies all the 1-specific parameters.
*
*  Alan Weinstein, ajw, 11/97
**********************************************************************

* Arguments:
      INTEGER IF1   ! input, flag for type of data required
      INTEGER IF2   ! input, flag for type of data required

* MC info
*#include "seq/clinc/qqpars.inc"
*#include "seq/clinc/qqprop.inc"
*#include "qqlib/seq/qqbrat.inc"

      INTEGER            JAK1,JAK2,JAKP,JAKM,KTOM
      COMMON / JAKI   /  JAK1,JAK2,JAKP,JAKM,KTOM
      REAL*4 RRR(1)
      REAL PARM(4,100)
      integer imixpp(300)
      INTEGER INIT,I,J
      REAL C1270,C1402,A1270_KSPI,A1270_KRHO,A1402_KSPI,A1402_KRHO
      REAL CG1,CG2,R,BRA1,BRKS
      SAVE INIT,PARM
      DATA INIT/0/

**********************************************************************
* Initialize return variable:
      PKORB = 0.

**********************************************************************
* Initialize:
      IF (INIT.EQ.0) THEN
        INIT = 1

C        CALL VZERO(PARM,400)
        DO I=1,4
        DO J=1,100
          PARM(I,J) = 0
        END DO
        END DO

C Youd better be using korb.dec, NOT decay.dec!!!!
C masses (needed in dist/inimas, formf/form*, etc)
        PARM(1, 1) = 1.777000   ! TAU
        PARM(1, 2) = 0.         ! NUTA
        PARM(1, 3) = 0.000511   ! EL
        PARM(1, 4) = 0.         ! NUEL
        PARM(1, 5) = 0.105658   ! MU
        PARM(1, 6) = 0.         ! NUMU
        PARM(1, 7) = 0.134976   ! PIZ
        PARM(1, 8) = 0.139570   ! PI+
        PARM(1, 9) = 0.769900   ! RHO+
        PARM(1,10) = 1.275000   ! A1+
        PARM(1,11) = 0.493677   ! K+
        PARM(1,12) = 0.497670   ! KZ
        PARM(1,13) = 0.891590   ! K*+
        PARM(1,14) = 0.781940   ! OMEG
        PARM(1,15) = 1.370000   ! RHOP+
        PARM(1,16) = 1.700000   ! K*P+
        PARM(1,17) = 1.461000   ! A1P+
        PARM(1,18) = 1.300000   ! PIP+
        PARM(1,19) = 1.270000   ! K1A+
        PARM(1,20) = 1.402000   ! K1B+
        PARM(1,21) = 1.465000   ! RHOPP+
        PARM(1,22) = 1.700000   ! RHOPPP+
	
C widths (needed in dist/inimas, formf/form*, etc)
        PARM(2, 1) = 0.          ! TAU
        PARM(2, 2) = 0.          ! NUTA
        PARM(2, 3) = 0.          ! EL
        PARM(2, 4) = 0.          ! NUEL
        PARM(2, 5) = 0.          ! MU
        PARM(2, 6) = 0.          ! NUMU
        PARM(2, 7) = 0.          ! PIZ
        PARM(2, 8) = 0.          ! PI+
        PARM(2, 9) = 0.1512      ! RHO+
        PARM(2,10) = 0.700       ! A1+
        PARM(2,11) = 0.          ! K+
        PARM(2,12) = 0.          ! KZ
        PARM(2,13) = 0.0498      ! K*+
        PARM(2,14) = 0.00843     ! OMEG
        PARM(2,15) = 0.510       ! RHOP+
        PARM(2,16) = 0.235       ! K*P+
        PARM(2,17) = 0.250       ! A1P+
        PARM(2,18) = 0.400       ! PIP+
        PARM(2,19) = 0.090       ! K1A+
        PARM(2,20) = 0.174       ! K1B+
        PARM(2,21) = 0.310       ! RHOPP+
        PARM(2,22) = 0.235       ! RHOPPP+

C Now store mixing parameters for 2pi and 4pi FFs 
C needed in tauola/fpik, tauola/bwigs, formf/form* , formf/curr :

        PARM(3,15) = -0.145

        IMIXPP(205)=1
        IMIXPP(207)=1
        IMIXPP(209)=1
        IMIXPP(211)=1
        IMIXPP(201)=1
        IMIXPP(203)=1
        IMIXPP(213)=1
        IMIXPP(215)=1


        IF (IMIXPP(205).NE.0) PARM(3,15) = -0.110
        IF (IMIXPP(207).NE.0) PARM(3,16) = -0.038
        IF (IMIXPP(209).NE.0) PARM(3,17) = 0.00
        IF (IMIXPP(211).NE.0) PARM(3,18) = 0.00
        IF (IMIXPP(201).NE.0) PARM(3,19) = 1.0
        IF (IMIXPP(203).NE.0) PARM(3,20) = 0.8
        IF (IMIXPP(213).NE.0) PARM(3,21) = -0.110
        IF (IMIXPP(215).NE.0) PARM(3,22) = -0.110

        PRINT *,' KORB: rho/rhop -> pi-pi0 mixing:'
        PRINT *,' KORB: rho   =',PARM(1,9) ,PARM(2,9)
        PRINT *,' KORB: rhop  =',PARM(1,15),PARM(2,15),PARM(3,15)
        PRINT *,' KORB: K*/K*prime -> Kpi mixing:'
        PRINT *,' KORB: kstp  =',PARM(1,16),PARM(2,16),PARM(3,16)
        PRINT *,' KORB: a1/a1prime -> 3pi, KKpi mixing:'
        PRINT *,' KORB: a1    =',PARM(1,10),PARM(2,10)
        PRINT *,' KORB: a1prim=',PARM(1,17),PARM(2,17),PARM(3,17)
        PRINT *,' KORB: K1A/K1B -> Kpipi mixing:'
        PRINT *,' KORB: K1A   =',PARM(1,19),PARM(2,19),PARM(3,19)
        PRINT *,' KORB: K1B   =',PARM(1,20),PARM(2,20),PARM(3,20)
        PRINT *,' KORB: rho/rhop/rhopp -> 4pi mixing:'
        PRINT *,' KORB: rho   =',PARM(1,9) ,PARM(2,9)
        PRINT *,' KORB: rhopp =',PARM(1,21),PARM(2,21),PARM(3,21)
        PRINT *,' KORB: rhoppp=',PARM(1,22),PARM(2,22),PARM(3,22)

C amplitudes for curr_cleo.F:
C for (3pi)-pi0: 4pi phase space; rho0pi-pi0; rho-pi+pi-; rho+pi-pi-; pi-omega
        PARM(3,31) = 0.
        PARM(3,32) = 0.1242
        PARM(3,33) = 0.1604
        PARM(3,34) = 0.2711
        PARM(3,35) = 0.4443
C for pi-3pi0: 4pi phase space; rho-pi0pi0
        PARM(3,36) = 0.
        PARM(3,37) = 1.0

C Modify amplitudes for 4pi form-factor in formf/curr, from korb.dec:
CCC        IF (IPLIST(2,282).EQ.5) THEN
        IPLIST=0
        IF (IPLIST.EQ.5) THEN
        PARM(3,31) = 0.0000
        PARM(3,32) = 0.1242
        PARM(3,33) = 0.1604
        PARM(3,34) = 0.2711
        PARM(3,35) = 0.4443
        PARM(3,36) = 0.0000
        PARM(3,37) = 1.0000
        END IF

        PRINT *,' KORB: 3PI-PI0 PARAMS:',(PARM(3,I),I=31,35)
        PRINT *,' KORB: PI-3PI0 PARAMS:',(PARM(3,I),I=36,37)

C The 4pi models are the most complicated in TAUOLA.
C If the user has not modified any parameters of the 4pi model,
C we can use the WTMAX determined with many trials.
        IF (ABS(PARM(3,31)-0.0000).GT.0.0001 .OR. 
     1      ABS(PARM(3,32)-0.1242).GT.0.0001 .OR.
     1      ABS(PARM(3,33)-0.1604).GT.0.0001 .OR.
     1      ABS(PARM(3,34)-0.2711).GT.0.0001 .OR.
     1      ABS(PARM(3,35)-0.4443).GT.0.0001 ) THEN
           PARM(3,38) = -1.
        ELSE
           PARM(3,38) = 6.9673671E-14
        END IF

        IF (ABS(PARM(3,36)-0.0000).GT.0.0001 .OR. 
     1      ABS(PARM(3,37)-1.0000).GT.0.0001 ) THEN
           PARM(3,39) = -1.
        ELSE
           PARM(3,39) = 3.5374880E-13
        END IF


C phases for curr_cleo.F:
        PARM(3,42) = -0.40
        PARM(3,43) =  0.00
        PARM(3,44) = -0.20+3.1416
        PARM(3,45) = -1.50

C rho' contributions to rho' -> pi-omega:
        PARM(3,51) = -0.10
        PARM(3,52) =  1.00
        PARM(3,53) = -0.10
        PARM(3,54) = -0.04

C rho' contribtions to rho' -> rhopipi:
        PARM(3,55) =  1.00
        PARM(3,56) =  0.14
        PARM(3,57) = -0.05
        PARM(3,58) = -0.05

C rho contributions to rhopipi, rho -> 2pi:
        PARM(3,59) =  1.000
        PARM(3,60) = -0.145
        PARM(3,61) =  0.000

C Set the BRs for (A1+ -> rho+ pi0) and (K*+ -> K0 pi+)
C needed in dist/taurdf:
        PARM(4,1) = 0.4920                 ! BRA1+
        PARM(4,2) = 0.4920                 ! BRA1-
        PARM(4,3) = 0.6660                 ! BRKS+
        PARM(4,4) = 0.6660                 ! BRKS-
        PARM(4,5) = 0.5                    ! BRK0
        PARM(4,6) = 0.5                    ! BRK0B

C amplitude coefficients for tau -> K1(1270) / K1(1402)
        C1270 = PARM(3,19)
        C1402 = PARM(3,20)
        IF (C1270.EQ.0.AND.C1402.EQ.0.) THEN
           C1270 = 1.
           C1402 = 0.6
        END IF
C From PDG96, square roots of branching fractions:
        A1270_KSPI = SQRT(0.16)
        A1270_KRHO = SQRT(0.42)
        A1402_KSPI = SQRT(0.94)
        A1402_KRHO = SQRT(0.03)
C C-G coefficients for K1- -> CG1 * |K- pi0> + CG2 * |K0bar pi->
        CG1 = -SQRT(2./3.)
        CG2 =  SQRT(1./3.)
C and the resulting amplitudes (times normalized FF):
        PARM(3,81) = C1270*A1270_KSPI*CG1  ! K1270 -> K*0B pi- 
        PARM(3,82) = C1402*A1402_KSPI*CG1  ! K1402 -> K*0B pi- 
        PARM(3,83) = C1270*A1270_KRHO*CG1  ! K1270 -> K0B rho- 
        PARM(3,84) = C1402*A1402_KRHO*CG1  ! K1402 -> K0B rho-
        PARM(3,85) = C1270*A1270_KSPI*CG2  ! K1270 -> K*- pi0
        PARM(3,86) = C1402*A1402_KSPI*CG2  ! K1402 -> K*- pi0
        PARM(3,87) = C1270*A1270_KRHO*CG2  ! K1270 -> K- rho0
        PARM(3,88) = C1402*A1402_KRHO*CG2  ! K1402 -> K- rho0

      END IF
**********************************************************************

      R = 0.
      IF (IF1.GE.1 .AND. IF1.LE.4 .AND. IF2.GE.1 .AND. IF2.LE.100) THEN
         R = PARM(IF1,IF2)

CAJW 4/4/94  Better to decide on A1 br now, avoid DADMAA/DPHSAA problem.
        IF (IF1.EQ.4.AND.JAK1.EQ.5) THEN 
          IF (IF2.EQ.11) THEN
C Return the BR used in the last call:
             R = BRA1
          ELSE IF (IF2.EQ.1) THEN
            BRA1 = R
            CALL RANMAR(RRR,1)
            IF (RRR(1).LT.BRA1) THEN
              R = 1.     ! 3pi
            ELSE 
              R = 0.     ! pi-2pi0
            END IF
            BRA1 = R
          END IF
        ELSEIF (IF1.EQ.4.AND.JAK1.EQ.7) THEN
          IF (IF2.EQ.13) THEN
C Return the BR used in the last call:
             R = BRKS
          ELSE IF (IF2.EQ.3) THEN
            BRKS = R
            CALL RANMAR(RRR,1)
            IF (RRR(1).LT.BRKS) THEN
              R = 1.     ! K0 pi-
            ELSE 
              R = 0.     ! K- pi0
            END IF
            BRKS = R
          END IF
        END IF

      END IF

      PKORB = R
      RETURN
      END
