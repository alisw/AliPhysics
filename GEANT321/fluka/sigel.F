*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:20:04  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.45  by  S.Giani
*-- Author :
*$ CREATE SIGEL.FOR
*COPY SIGEL
*
*=== sigel ============================================================*
*
      SUBROUTINE SIGEL ( IT, AA, EKIN, POO, SEL, ZL )
 
#include "geant321/dblprc.inc"
#include "geant321/dimpar.inc"
#include "geant321/iounit.inc"
C********************************************************************
C
C     Slightly changed on 11 february 1991  by  Alfredo Ferrari
C     (The kinetic energy has been added to the input variables,
C      a check is made on the minimum kinetic energy to have a
C      an elastic scattering with the present coding, sigmabar,
C      xsi, xsibar, omega, omegabar have been added with very
C      rough correspondence to the existing particles)
C
C     VERSION BY                     J. RANFT
C                                    LEIPZIG
C     LAST CHANGE 18. JULY 92 BY     A. FERRARI
C                                    INFN - MILAN
C
C
C     THIS IS A SUBROUTINE OF FLUKA82 TO GIVE ELASTIC SCATTERING
C     LENGTH AND CROSS SECTION
C
C     INPUT VARIABLES:
C     IT     = PARTICLE TYPE
C     AA     = ATOMIC WEIGHT OF THE NUCLEUS
C     EKIN   = PARTICLE KINETIC ENERGY
C     POO    = PARTICLE MOMENTUM IN GEV/C
C
C     OUTPUT VARIABLES:
C     SEL    = ELASTIC CROSS SECTION IN MB
C     ZL     = INTERACTION LENGTH IN G/CM**2
C
C     OTHER IMPORTANT VARIABLES:
C        SIG    = PROTON/NUCLEI CROSS SECTIONS
C        SEG    = PION/NUCLEI CROSS SECTIONS
C        P      = MOMENTUMS FOR WHICH THE CROSS SECTIONS ARE GIVEN IN
C                 SIG AND SEG
C        A      = NUCLEI FOR WHICH THE CROSS SECTIONS ARE GIVEN IN
C                 SIG AND SEG
C        PLAB   = MOMENTUMS FOR WHICH THE TOTAL CROSS SECTIONS ARE
C                 GIVEN IN SITO
C        SITO   = TOTAL HADRON NUCLEON CROSS SECTIONS FOR NUCLEONS,
C                 PIONS, KAONS AND ANTI-NUCLEONS.
C        ALP    =  EXPONENT OF THE PARAMETRIZATION FOR ANTI-PROTONS,
C                  RANTI-NEUTRONS AND KAONS
C        BET    =  MULTIPLIER OF PARAMETRIZATION FOR ANTI-PROTONS,
C                  ANTI-NEUTRONS AND KAONS
C
C     NOTE1: PRESENTLY CROSS SECTIONS ARE ASSUMED TO BE CONSTANT
C     ABOVE 10.0 GEV/C FOR ALL PARTICLES AND
C     BELOW 0.3 GEV/C FOR NUCLEONS AND BELOW 0.13 GEV/C FOR PIONS
C
C     NOTE2: FOR HADRONS OTHER THAN (1=PROTON,2=ANTI PROTON,8=
C     NEUTRON,9=ANTI NEUTRON,13=POSITIVE PION,14=NEGATIVE PION,15=
C     POSITIVE KAON,16=NEGATIVE KAON,24=NEUTRAL KAON,25=NEUTRAL ANTI
C     KAON) SEE TABLE ITT TO SEE THE CORRESPONDANCE
C
C     NOTE3: FOR LEPTONS AND PHOTONS PRACTICALLY ZERO CROSS SECTION
C     IS RETURNED.
C
C********************************************************************
C
#include "geant321/paprop.inc"
      PARAMETER (AVOGMB=AVOGAD*1.D-27)
C--------------------------------------------------------------------
      DIMENSION SIG(13,9),SEG(16,9),P(16),A(9),ITT(39)
      DIMENSION PLAB(19),SITO(19,4),ALP(3),BET(3)
      DIMENSION REA(9,9),STOT(9)
      SAVE A, P, SIG, SEG, ITT, SITO, PLAB, ALP, BET, STOT, REA
      DATA A/9.D0,12.D0,27.D0,47.9D0,55.9D0,63.5D0,112.4D0,
     &207.2D0,238.1D0/
      DATA P/.13D0,.19D0,.25D0,.3D0,.4D0,.5D0,.6D0,.8D0,1.D0,
     &1.5D0,2.D0,3.D0,4.D0,5.D0,6.D0,10.D0/
      DATA SIG/ 485.D0,223.D0,112.D0,82.D0,66.D0,78.D0,96.D0,102.D0,
     &100.D0,98.D0,95.D0,90.D0,79.D0,
     (680.D0,348.D0,175.D0,103.D0,84.D0,87.D0,106.D0,112.D0,111.D0,
     &108.D0,107.D0,105.D0,101.D0,
     (1200.D0,738.D0,387.D0,196.D0,191.D0,200.D0,248.D0,264.D0,
     &264.D0,257.D0,252.D0,247.D0,228.D0,
     (1658.D0,1110.D0,635.D0,364.D0,332.D0,356.D0,404.D0,408.D0,
     &407.D0,404.D0,398.D0,396.D0,384.D0,
     (1730.D0,1270.D0,725.D0,400.D0,375.D0,412.D0,495.D0,505.D0,
     &495.D0,492.D0,487.D0,485.D0,475.D0,
     (1875.D0,1470.D0,835.D0,480.D0,450.D0,450.D0,535.D0,
     &580.D0,555.D0,540.D0,535.D0,530.D0,
     (525.D0,2040.D0,2160.D0,1335.D0,850.D0,740.D0,760.D0,880.D0,
     &905.D0,860.D0,840.D0,820.D0,815.D0,800.D0,
     (2340.D0,2980.D0,2270.D0,1450.D0,1230.D0,1230.D0,
     &1380.D0,1420.D0,1410.D0,1380.D0,1360.D0,
     (1350.D0,1320.D0,2680.D0,3220.D0,2530.D0,1630.D0,1420.D0,1450.D0,
     &1570.D0,1600.D0,1590.D0,1575.D0,1560.D0,1550.D0,1540.D0/
      DATA SEG/24.D0,128.D0,249.D0,256.D0,202.D0,124.D0,73.D0,60.D0,
     &64.D0,69.D0,62.D0,50.D0,44.D0,42.D0,42.D0,41.D0,21.D0,156.D0,
     (273.D0,280.D0,220.D0,212.D0,94.D0,80.D0,82.D0,85.D0,80.D0,73.D0,
     (69.D0,67.D0,66.D0,64.D0,56.D0,296.D0,560.D0,574.D0,467.D0,350.D0,
     &235.D0,210.D0,210.D0,200.D0,190.D0,183.D0,176.D0,170.D0,165.D0,
     (155.D0,100.D0,500.D0,895.D0,880.D0,690.D0,520.D0,378.D0,
     (355.D0,384.D0,373.D0,352.D0,320.D0,300.D0,288.D0,280.D0,
     &262.D0,75.D0,500.D0,965.D0,990.D0,775.D0,525.D0,410.D0,410.D0,
     (433.D0,440.D0,425.D0,395.D0,374.D0,355.D0,340.D0,303.D0,125.D0,
     (570.D0,1025.D0,1100.D0,825.D0,575.D0,418.D0,458.D0,500.D0,
     &480.D0,460.D0,440.D0,422.D0,400.D0,384.D0,355.D0,300.D0,880.D0,
     (1480.D0,1550.D0,1380.D0,940.D0,710.D0,720.D0,810.D0,760.D0,
     (740.D0,700.D0,665.D0,645.D0,620.D0,570.D0,550.D0,1475.D0,
     &2250.D0,2350.D0,1850.D0,1500.D0,
     (1120.D0,1210.D0,1480.D0,1440.D0,1400.D0,1320.D0,
     &1250.D0,1210.D0,1170.D0,1065.D0,540.D0,
     (1300.D0,2220.D0,2560.D0,1980.D0,1650.D0,1160.D0,
     &1360.D0,1600.D0,1560.D0,1510.D0,1410.D0,
     (1350.D0,1300.D0,1270.D0,1200.D0/
C     DATA ITT/1,7,0,0,0,0,0,2,8,0,0,9,3,4,6,5,1,2,9,1,1,1,3,9,10,
      DATA ITT/1,7,0,0,0,0,0,2,8,0,0,9,3,4,6,5,2,8,9,1,1,2,3,9,10,
     &         3,0,0,0,0,7,2,7,2,8,1,7,1,7/
      DATA PLAB/.3D0,.4D0,.5D0,.6D0,.7D0,.8D0,.9D0,1.D0,1.1D0,
     &1.2D0,1.3D0,1.4D0,1.5D0,2.D0,3.D0,4.D0,
     *5.D0,6.D0,10.D0/
      DATA SITO/66.8D0,63.6D0,40.35D0,31.25D0,31.1D0,
     *35.1D0,36.7D0,44.15D0,38.3D0,33.25D0,
     *29.75D0,29.3D0,29.95D0,26.55D0,24.6D0,22.95D0,
     *22.75D0,22.95D0,21.55D0,
     *12.5D0,14.1D0,13.5D0,12.75D0,12.85D0,13.9D0,15.6D0,
     *17.25D0,18.9D0,19.5D0,18.95D0,18.85D0,
     *18.45D0,18.2D0,17.5D0,17.7D0,17.5D0,17.25D0,17.4D0,
     *39.65D0,38.75D0,26.9D0,22.D0,22.D0,24.5D0,26.15D0,30.7D0,28.6D0,
     &26.4D0,24.35D0,24.1D0,24.2D0
     (,22.4D0,21.05D0,20.3D0,20.1D0,20.1D0,19.5D0,
     (280.D0,199.7D0,171.1D0,154.3D0,140.D0,130.D0,116.8D0,117.4D0,
     &111.6D0,109.D0,106.5D0,
     (102.8D0,100.D0,90.2D0,76.7D0,68.D0,62.8D0,60.7D0,56.D0/
      DATA ALP/0.823D0,0.843D0,0.630D0/
      DATA BET/1.26D0,1.31D0,0.90D0/
      DATA STOT /15.D0,20.D0,30.D0,40.D0,60.D0,80.D0,
     &100.D0,150.D0,200.D0/
      DATA REA / .20D0,.23D0,.27D0,.30D0,.35D0,.40D0,.47D0,.55D0,.60D0,
     2           .22D0,.26D0,.31D0,.35D0,.40D0,.45D0,.51D0,.59D0,.63D0,
     3           .24D0,.29D0,.36D0,.42D0,.50D0,.56D0,.60D0,.66D0,.68D0,
     4           .26D0,.32D0,.42D0,.49D0,.58D0,.63D0,.66D0,.71D0,.72D0,
     5           .27D0,.33D0,.44D0,.51D0,.61D0,.65D0,.68D0,.72D0,.74D0,
     6           .28D0,.35D0,.46D0,.53D0,.63D0,.66D0,.69D0,.73D0,.745D0,
     7           .35D0,.42D0,.53D0,.62D0,.69D0,.72D0,.74D0,.77D0,.78D0,
     8           .42D0,.51D0,.62D0,.69D0,.75D0,.77D0,.79D0,.81D0,.82D0,
     9           .44D0,.53D0,.64D0,.70D0,.76D0,.78D0,.80D0,.81D0,.82D0 /
C
C
C
      SEL = AZRZRZ
      ZL  = AINFNT
      IF(AA.LT.0.99D0)RETURN
      IPOL=0
      PO=POO
      EKE=EKIN
* Set a very large kinetic energy for the call to Nizl to bypass
* any check on inelastic events thresholds
      EKINFN=AINFNT
      IIT=ITT(IT)
      IF (IIT.EQ.0) RETURN
C---------------------------------------------------------
C**                          ELASTIC SCATTERING ON PROTONS
C         HJM 10/88             REASONABLE FOR P, N, PI+/-
C**
      IF((AA.LT.1.5D0).AND. (IT.EQ.1.OR.IT.EQ.8.OR.IT.EQ.13
     *  .OR.IT.EQ.14.OR.IT.EQ.2.OR.IT.EQ.9)) THEN
C**      EKE=SQRT(PO**2+AM(IT)**2) - AM(IT)
* Return if kinetic energy is below 15 MeV (A. Ferrari)
         IF ( EKE .LT. 0.015D+00 ) RETURN
         IF ( IT .EQ. 9 ) THEN
            IJT = 1
         ELSE IF ( IT .EQ. 2 ) THEN
            IJT = 8
         ELSE
            IJT = IT
         END IF
         CALL SIHAEL(IJT,EKE,PO,AA,SEL)
         GOTO 124
      END IF
C**
C                            NEUTRON-NUCLEUS ELASTIC SCATTERING
C                            DATA FROM HETKFA2 FOR  EKIN .GT. 15 MEV
C                            FOR PLOTS SEE
C                               P. CLOTH ET AL.,
C                               HERMES - A MC PROGRAM SYSTEM ...
C                               JUEL-2203 (MAY 1988)
C
*     IF((IT.EQ.8.OR.IT.EQ.1.OR.IT.EQ.2.OR.IT.EQ.9).AND.PO.LT.20.D0)
*    &   THEN
      IF(IT.EQ.8.AND.PO.LT.20.D0)THEN
* Return if kinetic energy is below 15 MeV (A. Ferrari)
         IF ( EKE .LT. 0.015D+00 ) RETURN
         IF ( IT .EQ. 9 ) THEN
            IJT = 1
         ELSE IF ( IT .EQ. 2 ) THEN
            IJT = 8
         ELSE
            IJT = IT
         END IF
         IF(PO.GT.10.D0) THEN
            IPOL=1
            PO=10.D0
            EKE=SQRT(PO**2+AM(IT)**2) - AM(IT)
         END IF
C**      EKE=SQRT(PO**2+AM(IT)**2) - AM(IT)
         CALL SIHAEL(IJT,EKE,PO,AA,SEL)
         IF(IPOL.EQ.1) GOTO 210
         GOTO 124
      ENDIF
      IF ( EKE .LT. 0.020D+00 ) RETURN
C-----------------------------------------------------------
C
C
C********************************************************************
C     CALCULATE THE NEW PARTICLE NUMBER IIT:   1=P,2=N,3=PI+,4=PI-,
C     5=K-,6=K+,7=P BAR,8=N BAR,9=K ZERO ,10=K ZERO BAR
C********************************************************************
C
      IF((IIT.EQ.7).OR.(IIT.EQ.8)) GOTO 200
      IF (PO.GT.20.D0) GOTO 200
      IF(PO.LE.10.D0) GOTO 30
      PO=10.D0
      IPOL=1
 30   CONTINUE
      IF(IIT.LE.4) GO TO 10
C
C********************************************************************
C     MOMENTUM INDEX K FOR KAONS ANTI KAONS AND ANTI NUCLEONS
C********************************************************************
C
      DO 3 K=1,19
      IF(PO.LE.PLAB(K)) GO TO 40
    3 CONTINUE
      K=19
   40 GO TO 7
C
C********************************************************************
C     CALCULATE THE MOMENTUM INDEX K FOR NUCLEONS AND PIONS
C     CALCULATE THE MASS INDEX J OF THE NUCLEUS
C********************************************************************
C
   10 CONTINUE
      DO 22 K=1,16
      IF(PO.LE.P(K)) GO TO 23
   22 CONTINUE
      K=16
   23 CONTINUE
      DO 5 I=2,8
      IF(AA.LE.A(I)) GO TO 6
      GO TO 5
    6 CONTINUE
      J=I-1
      GO TO 7
    5 CONTINUE
      J=8
C
C********************************************************************
C     SELECT THE FORMULEI TO BE USED FOR DIFFERENT PARTICLE TYPES
C********************************************************************
C
    7 CONTINUE
C             P , N ,PI+,PI-,K- ,K+ ,AP ,AN ,K0 ,AK0
      GO TO (101,101,113,113,116,115,102,109,115,116),IIT
C******************** PROTONS,NEUTRONS,OTHERS
  101 K=K-3
      IF(K.LT.1) K=1
      ALOGA=LOG(A(J+1)/A(J))
      AAA=AA/A(J)
      SI1=SIG(K,J)* AAA     **(LOG(SIG(K,J+1)/SIG(K,J))/ALOGA)
      IF(K.EQ.1) GO TO 2000
      KK=K-1
      SI2=SIG(KK,J)* AAA     **(LOG(SIG(KK,J+1)/SIG(KK,J))/ALOGA)
      K=K+3
      KK=KK+3
      SI=SI1+(PO-P(K))*(SI2-SI1)/(P(KK)-P(K))
C
      SEL=SI
C
      GO TO 121
C******************** CHARGED PIONS
  113 CONTINUE
      ALOGA=LOG(A(J+1)/A(J))
      AAA=AA/A(J)
      SI1=SEG(K,J)* AAA     **(LOG(SEG(K,J+1)/SEG(K,J))/ALOGA)
      IF(K.EQ.1) GO TO 2000
      KK=K-1
      SI2=SEG(KK,J)* AAA     **(LOG(SEG(KK,J+1)/SEG(KK,J))/ALOGA)
      SI=SI1+(PO-P(K))*(SI2-SI1)/(P(KK)-P(K))
C
      SEL=SI
C
      GO TO 121
C******************** K-,K ZERO BAR
  116 CONTINUE
      IA=1
      IS=1
      GO TO 122
C******************** K+, K ZERO
  115 CONTINUE
      IA=2
      IS=2
      GO TO 122
C******************** P BAR
  102 CONTINUE
C******************** N BAR
  109 CONTINUE
C
      PO=POO
      GOTO 200
C
C
C********************************************************************
C     KAONS, ANTI KAONS
C********************************************************************
C
  122 KK=K-1
      IF(K.EQ.1) GO TO 140
      PKK=PLAB(KK)
      SIKK=SITO(KK,IS)
      SI=(SITO(K,IS)-SIKK)*(PO-PKK)/(PLAB(K)-PKK)+SIKK
      GO TO 141
  140 SI=SITO(K,IS)
  141 SI1=SI
      SI=BET(IA)*SI*AA**ALP(IA)
      IV=IT
      IF(IV.NE.24)GO TO 150
      IV=15
      SI=SI*2.06D0
      GO TO 151
 150  IF(IV.EQ.25) IV=16
 151  CALL NIZL(IV,AA,EKINFN,PO,SINEL,ZLIN)
C
      SEL=SI-SINEL
      IF(IPOL.EQ.1) GOTO 210
      GOTO 121
C
C********************************************************************
C     AND NOW THE SCATTERING LENGTH IN G/CM**2
C********************************************************************
C
 121  CONTINUE
      IF(IPOL.EQ.1) GOTO 210
 124  CONTINUE
C
      IF (SEL.LT.ANGLGB) THEN
         SEL = AZRZRZ
         ZL  = AINFNT
      ELSE
         ZL=AA/(AVOGMB*SEL)
      END IF
      RETURN
C
C********************************************************************
C     WE ARE IN THE LOWEST MOMENTUM BIN
C********************************************************************
C
 2000 CONTINUE
      SI=SI1
      SEL=SI
      GO TO 121
C***
C   ENTRY FOR SMOOTHING OF SIGEL BETWEEN 10. AND 20. GEV/C
C***
 210  CONTINUE
      PO=20.D0
C***
C   APPROXIMATION FOR HIGH ENERGIES
C***
 200  CONTINUE
C
      IT1=IT
      IF((IT.EQ.2).OR.(IT.EQ.9)) IT1=1
C
      STO=SHPTOT(IT1,PO)
C
C   MASS NUMBER INDEX
C***
      DO 201 IA=2,8
      IF(AA.GT.A(IA)) GOTO 201
      JA=IA-1
      GOTO 202
 201  CONTINUE
      JA=8
 202  CONTINUE
C***
C   SIGTOT INDEX
C***
      DO 203 IS=2,8
      IF(STO.GT.STOT(IS)) GOTO 203
      JS=IS-1
      GOTO 204
 203  CONTINUE
      JS=8
 204  CONTINUE
C
      DA1=A(JA+1)-A(JA)
      DA2=AA - A(JA)
      RR=REA(JS,JA)
      R1=RR + DA2*(REA(JS,JA+1)-RR)/DA1
      RR=REA(JS+1,JA)
      R2=RR + DA2*(REA(JS+1,JA+1)-RR)/DA1
      RACT=R1 + (STO-STOT(JS))*(R2-R1)/(STOT(JS+1)-STOT(JS))
C
      CALL NIZL(IT,AA,EKINFN,PO,SINEL1,ZLIN)
      SEL1=RACT*SINEL1
      IF(IPOL.EQ.1) GOTO 211
      SEL=SEL1
      GOTO 124
C
 211  CONTINUE
      SEL=SEL + (SEL1-SEL)*(POO-10.D0)/10.D0
      GOTO 124
C
C
C********************************************************************
C     FORMATS
C********************************************************************
C
 1000 FORMAT('          WARNING AT CALL SIGEL  ',I5)
      END
