*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:19:57  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.43  by  S.Giani
*-- Author :
*$ CREATE NUCRIV.FOR
*COPY NUCRIV
*
*=== nucriv ===========================================================*
*
      SUBROUTINE NUCRIV (ITTTT,ELAB,CX,CY,CZ,BBTAR,ZZTAR,RHOO)
 
#include "geant321/dblprc.inc"
#include "geant321/dimpar.inc"
#include "geant321/iounit.inc"
*
*----------------------------------------------------------------------*
*                                                                      *
*  Last change  on  16-sep-92    by    Alfredo Ferrari, Infn - Milan   *
*                                                                      *
*  Nucrin90 by A. Ferrari: tentative version to overcome many troubles *
*                          in energy, momentum and charges conservation*
*                                                                      *
*----------------------------------------------------------------------*
*
C
C--------------------------------------------------
C*** FINUC: NUCRIN FINAL STATE PARTICLE LIST WITH KINEM. VARIABLES
C*** FINLSP HADRIN FINAL STATE PARTICLE LIST WITH KINEM. VARIABLES
C*** (NUMBER OF PARTICLES,PARTICLE TYPE INDEX,DIRECTION COSINES,ENERGY
C*** ABSOLUTE MOMENTUM AND IN FINUC IN ADDIT. EXCITATION ENERGY
C--------------------------------------------------
#include "geant321/balanc.inc"
#include "geant321/corinc.inc"
#include "geant321/finuc2.inc"
#include "geant321/finlsp2.inc"
#include "geant321/nucdat.inc"
#include "geant321/resnuc.inc"
*     PARAMETER ( TVEPSI = 1.D-05 )
      COMMON / FKNUCF / DELEFT, EKRECL, V0EXTR, ITTA, ITJ, LVMASS
      DIMENSION EXSOP (2)
      REAL RNDM(1),RANGS1,RANGS2
      COMMON /FKFERM/ELABKE,A,TPNVK,ELATE,IFERT
C
C--------------------------------------------------
C***
C***
C***  SAMPLING OF A HADRON-NUCLEUS COLLISION EVENT
C***
C   ITTTT - TYPE OF INCOMING HADRON
C
C--------------------------------------------------
C***  POSSIBLE: N, AN,PI,K,AK,Y,AY BY THE INDEX ITTTT=IT=1,2,8-25
C   PLAB - MOMENTUM OF INCOMING HADRON
C
C***  RANGE UP TO 5. GEV/C (FROM ABOUT 0. ... 0.1 GEV/C)
C   ELAB - TOTAL ENERGY OF INCOMING HADRON IN GEV
C   CX,CY,CZ - DIRECTION COSINES
C***  ANUC,ZNUC = NUMBERS OF NUCLEONS AND PROTONS
C*** RHOO = MATERIAL DENSITY (G/CM**3)
C
C*** FINAL STATE PARTICLE CHARACTERISTICS TABLE IN /FINUC/:
C*** IRN   - NUMBER OF FINAL STATE PARTICLES
C*** ITRN(I) - FINAL STATE PARTICLE TYPE INDEX
C*** CXRN,CYRN,CZRN (I) - DIRECTION COSINES OF F.S.P. (LAB SYST.)
C*** ELR,PLR (I) - LAB.ENERGY AND MOMENTUM OF F.S.P. (GEV, GEV/C)
C*** TV    - EXCITATION ENERGY (GEV)
C
C--------------------------------------------------
      DIMENSION THRESR(30)
      SAVE THRESR,ITPRF,IAMC,JAMC,INS,IXPI,JNUC
      DATA THRESR/1.9D0, 0.D0, 5*9.D0, 1.9D0, 0.D0, 3*9.D0, 1.08D0,
     * 1.08D0, 1.44D0, 1.08D0, 6*9.D0, 1.08D0,
     * 1.44D0, 1.08D0, 5*9.D0/
      DIMENSION IXPI (30)
*
*----------------------------------------------------------------------*
*                                                                      *
*     Ixpi :                                                           *
*             = -1 for antinucleons ( pbar, nbar, and for Lambdabar,   *
*                  see below )                                         *
*             =  0 for pi-, K- ( stopping particles )                  *
*             =  1 for other mesons and particles (p, n, pi+, pi0, K+, *
*                  K0, K0bar )                                         *
*             =  2 for leptons (e-, e+, nu, nubar, mu+, mu-, photons,  *
*                  Klong, Kshort )                                     *
*             =  4 for hyperons ( Lambda, Lambdabar, Sigma-, Sigma+,   *
*                  Sigma0 ). Actually for Lambdabar the inxpi flag in  *
*                  the code, which corresponds to the ixpi one of the  *
*                  projectile, is then set to - 1 as for antinucleons  *
*                                                                      *
*----------------------------------------------------------------------*
*
      DATA IXPI/1,-1,5*2,1,-1,3*2,1,0,1,0,4,4,2,3*4,1,1,1,5*2/
      COMMON /FKPERC/ IPERCO
      COMMON /FKENCO/ ETEST,TNKTE
C
C--------------------------------------------------
C***  PARTICLE CHARACTERISTICS: MASSES, DECAY WIDTH, LIFE TIME,ELECT.AND
C***  BARYONIC CHARGE, DECAY CHANNEL INDICEES
C--------------------------------------------------
*
*  /Abltis/ common is initialized in Hadden. Masses are the same as in
*  /Part/ common at least up to particle n. 94, the same for the baryonic
*  charge and the electric charge : what about particle 26????????????
*
      COMMON / FKABLT / AM (110), GA (110), TAU (110), ICH (110),
     &                  IBAR (110), K1(110), K2(110)
      COMMON / FKNUCT / ETHR, PTHR
      DIMENSION INS(30)
      DIMENSION AMHH(15),IAMC(15),JAMC(30)
      DATA IAMC/1,2,8,9,12,15,16,17,18,19,20,21,22,24,25/
      DATA JAMC/2*1,5*0,2*1,2*0,1,2*0,16*0/
      DATA INS/13,13,5*32,14,14,3*32,10,10
     *,12,11,7*32,15,15,5*32/
      DATA JNUC/0/
      DIMENSION ITPRF(110)
      DATA ITPRF/-1,-1,5*1,-1,-1,1,1,1,-1,-1,-1,-1,6*1,-1,-1,-1,85*1/
      LOGICAL LEMINU, LPRONU, LPCASC, LNCASC, LCORIN, LDELTX, LVMASS,
     &        LHYPER
*
      LCORIN = .FALSE.
      LVMASS = .FALSE.
      ANUC = BBTAR
      ZNUC = ZZTAR
      TPNVK  = 0.D+00
      IRNTOT = 0
      INNUC=1
      ITNUCR=ITTTT
      SICO   = 1.D0
      ELPH   = 0.D0
      PLABCO = 1.D0
C
C--------------------------------------------------
C*** SICO,PLABCO - EFFECTIVE CROSS SECTION- AND EFFECT. LABMOMENTUM
C***      CORRECTION FACTORS FOR Y-A- AND AY-A-COLLISIONS
C*** ITNUCR - HYPERON NUCLEUS COLL. PARTICLE TYPE FOR USE IN NIZL,
C***      IT - USED IN HADRIN
C--------------------------------------------------
      IBAT = 0
C--------------------------------------------------
C*** INXPI:
C*** ORDERING INDICEES FOR STOPPING PARTICLES (0), OTHER MESONS (1),
C*** ANTINUCLEONS (-1), HYPERONS (4),LEPTONS,KO,AKO(2),OTHER PART. (1)
C
C--------------------------------------------------
C*** INCLUSION OF HYPERON-NUCLEUS-COLLISIONS
C*** PRELIMINARY FOR SIGMA+ NO STRANGENES-CONSERVATION
C
C--------------------------------------------------
*
*  From here it is the number of the incoming nucleon
*
      IT = ITTTT
*
*  Nxpi = 0 means hadrin call is possible
*
      NXPI  = 0
      INXPI = IXPI(IT)
      IF ( IT .EQ. 23 .AND. ELAB - AM(IT) .LE. 0.05D+00 ) INXPI = 0
*  +-------------------------------------------------------------------*
*  |
      IF ( INXPI .EQ. 4 ) THEN
         CALL HYPERO ( IT, ITNUCR, SICO, PLABCO )
         LHYPER = .TRUE.
*  |
*  +-------------------------------------------------------------------*
*  |
      ELSE
         LHYPER = .FALSE.
      END IF
*  |
*  +-------------------------------------------------------------------*
C
C--------------------------------------------------
C*** CUT OFF ENERGY CONSTANTS:
C--------------------------------------------------
      ETHR  = 0.001D+00
C
C--------------------------------------------------
C*** CDDT,CEET PARAMETERS FOR GAUSSIAN WIDTH'
C--------------------------------------------------
      N = ITTTT
*
*  Particle 18 is Alambda
*
      IF ( ITTTT .EQ. 18 ) INXPI = -1
C
C--------------------------------------------------
C*** ETEST = ENERGY CONSERVATION TEST VARIABLE (SHOULD BE 0 AT RETURN)
C*** TLAB = KINETIC ENERGY
C--------------------------------------------------
      JJ    = JNUC
      CDDT  = 0.5D0
      IHACA = 0
      IN    = INS(N)
C
C--------------------------------------------------
C*** JJ=1:
C*** OPTION 1: GO TO NIZL, CROSS SECTION CALCULATION
C*** JJ=2:
C*** OPTION 2: EVENT SAMPLING
C--------------------------------------------------
      TLAB  = ELAB - AM (IT)
*  +-------------------------------------------------------------------*
*  |  Decide which is the target nucleon ( flag Itta ) : a proton .....
      CALL GRNDM(RNDM,1)
      IF ( RNDM (1) .LT. ZNUC / ANUC ) THEN
         ITTA = 1
         ITJ  = 1
*  |
*  +-------------------------------------------------------------------*
*  |                                            .... or a neutron
      ELSE
         ITTA = 8
         ITJ  = 2
      END IF
*  |
*  +-------------------------------------------------------------------*
      V0OLD = V0WELL (ITJ)
      PLAB = SQRT ( TLAB * ( ELAB + AM (IT) ) )
*  +-------------------------------------------------------------------*
*  |
      IF ( ABS ( PLAB - 5.D0 ) .GE. 4.99999D0 ) THEN
         WRITE(LUNOUT,99996) PLAB
         WRITE(LUNERR,99996) PLAB
*        STOP         Commented out A. Fasso' 1989
         IRN = NP0
         RETURN
99996    FORMAT (3(5H ****/),' Nucrin: projectile momentum',
     &          ' outside of the allowed region, plab, kproj:',
     *   1P,E15.5,3X,0P,I3,/,3(5H ****/))
      END IF
*  |
*  +-------------------------------------------------------------------*
*     IF ( JJ .EQ. 1 ) GO TO 1000
*-->-->-->-->--> go to cross section calculation if this option was se-
*                lected
C
C--------------------------------------------------
C*** IF LOW KINETIC ENERGY ( TLAB < ETHR ): ONLY EXCITATION
C*** STOPPING PARTICLES
C*** AND ANNIHILATION , INXPI=IXPI(IT) LE 0
C--------------------------------------------------
*  +-------------------------------------------------------------------*
*  |  Inxpi > 0: p, n, pi+, K+, Lambda, Sigma+, Sigma-, Sigma0
*  |  ( no lepton or pi0 projectile of course, now also pi0
*  |    are possible projectiles )
      IF (INXPI .GT. 0) THEN
C
C--------------------------------------------------
C*** ETHRR=THRESHOLD VALUE, FOR TLAB<ETHRR IS NO H-N COLLISION POSSIBLE,
C*** ONLY EXCITATION + CASCADE AND/OR ANNIHILATION CORRESP TO THE HADRN
*  |  +----------------------------------------------------------------*
*  |  |  Raised to 10 MeV: A. Ferrari, 11-apr-1991
         IF (TLAB .LE. 10.D+00 * ETHR) THEN
*  |  |   Here particles not stopping and not to be annihilated below
*  |  |   threshold
            ANOW  = ANOW  + IBAR (IT)
            KTARN = KTARN - IBAR (IT) + ICH (IT)
            ZNOW  = ZNOW  + ICH  (IT)
            KTARP = KTARP - ICH  (IT)
            AMNRES = AMUAMU * ANOW + 1.D-03 * FKENER ( ANOW, ZNOW ) -
     &               ZNOW * AMELEC + ELBNDE ( NINT (ZNOW) )
            TV  = ETTOT - AMNRES
            IRN = NP0
            RETURN
         END IF
*  |  |
*  |  +----------------------------------------------------------------*
*  |
*  +-------------------------------------------------------------------*
*  |    Here for ap, an, pi-, k-, alambda
      ELSE
*  |  +----------------------------------------------------------------*
*  |  |  Threshold is raised for particles with ixpi < 0 (ap, an,
*  |  |  alambda) to 80 MeV
         IF (INXPI .LT. 0) THEN
            ETHRR = 0.08D+00
*  |  |
*  |  +----------------------------------------------------------------*
*  |  |  for pi-, K- the threshold is Ethr = 1 MeV
         ELSE
            ETHRR = ETHR
         END IF
*  |  |
*  |  +----------------------------------------------------------------*
C
C--------------------------------------------------
C*** UP TO ETHRR=.08 GEV NO AN, AP IN FINAL STATE
C--------------------------------------------------
         AIO  = 10.D0
         XPI  = 1.D0 - AIO * ( TLAB - ETHRR )
*  |  Patch for pseudo pi0
         IF ( IT .EQ. 23 ) XPI = 1.D+00
*  |  For pi- and K- :
*  |                   Xpi < 0 for Tlab > 1/Aio + Ethrr = 101 MeV
*  |                   Xpi > 1 for Tlab < Ethrr = 1 MeV
*  |  For antibaryons:
*  |                   Xpi < 0 for Tlab > 1/Aio + Ethrr = 180 MeV
*  |                   Xpi > 1 for Tlab < Ethrr = 80 MeV
*  |  +----------------------------------------------------------------*
*  |  |   Nxpi=0 means hadrin-call is possible, else impossible
*  |  |   So Hadrin call is always possible for pi- and K- if
*  |  |   Tlab > 101 MeV, and for antibaryons if Tlab > 180 MeV.
*  |  |   Hadrin call is always impossible for pi- and K- if
*  |  |   Tlab < 1 MeV (never possible) and for antibaryons if
*  |  |   Tlab < 80 MeV. Actually Nxpi > 1 for antibaryons means forced
*  |  |   annihilation, but Hadrin can be called the same, even though
*  |  |   only final states with annihilation will be accepted
         CALL GRNDM(RNDM,1)
         IF ( RNDM(1) .GT. XPI ) THEN
            NXPI = 0
*  |  |                       particle non stopping (annihilating)
*  |  +----------------------------------------------------------------*
*  |  |                       particle stopping (annihilating)
         ELSE
            NXPI = 1
         END IF
*  |  |
*  |  +----------------------------------------------------------------*
      END IF
*  |
*  +-------------------------------------------------------------------*
*  +-------------------------------------------------------------------*
*  |    Incident pi- or K- stopping in the nucleus, no possibility
*  |    to call hadrin: it must be updated using the new stopping
*  |    particle module as soon as it will be available!!
      IF ( NXPI .GT. 0 .AND. IBAR (IT) .EQ. 0 ) THEN
*  |  Give the total lab energy of stopping pi-, K- into cascade and
*  |  exit: label 405 is where this is performed.
*  |  Before change a proton into a neutron to account for charge
*  |  conservation and adjust the total energy with the mass
*  |  difference (Note that for K- no strangeness conservation occurs
         IF ( IT .NE. 23 ) THEN
            ZNOW  = ZNOW - 1.D+00
            KTARP = KTARP + 1
            KTARN = KTARN - 1
         END IF
         TPK  = 0.D0
         TNK  = 0.D0
         TV   = 0.D0
         TOEFF = ELAB
         LCORIN = .TRUE.
         CALL CORSTP ( TOEFF )
         CALL GRANOR(RANGS1, RANGS2)
         TNKTE = MAX ( 0.003D+00, EKMNNU (2)
     &         * ( 0.5D+00 - 0.25D+00 * RANGS1 ) )
         TPKTE = MAX ( 0.003D+00, EKMNNU (1)
     &         * ( 0.5D+00 - 0.25D+00 * RANGS2 ) )
         GO TO 405
*  |
*  +-->-->-->-->--> go to the cascade simulation
      END IF
*  |
*  +----------------------------------------------------------------*
C
C--------------------------------------------------
C*** ELABKE = INVARIANT KINETIC ENERGY OF THE H-A-SYSTEM+PROJ.H.MASS
C*** A = TOTAL ENERGY OF THE PRIMARY PARTICLE IN LABSYSTEM
C*** IFORBI,IFERT LOOP COUNTING INDICEES FOR CHOICE OF FERMION MOMENTUM
C*** AMNTAR - NUCLEUS MASS
C*** UMOJAN - H-A-C.M.S.-ENERGY
C*** ELABKE - KINETIC H-A-ENERGY IN C.M.S
C--------------------------------------------------
*
      UMOJAN = SQRT ( AM(IT)**2 + AMNTAR**2 + 2.D+00 * ELAB * AMNTAR )
C
C--------------------------------------------------
C***  START WITH THE EVENT SIMULATION
C--------------------------------------------------
      ELABKE = UMOJAN - AMNTAR
      A = ELAB
      IFORBI = 0
      IFERT  = 0
*  Call Ferset for setting up anything for the Fermi motion
      CALL FERSET
      AMMIT = AM (IT)
      AMMTA = AM (ITTA)
      AIT2  = AMMIT**2
      ATA2  = AMMTA**2
*  Ecmo is the square of the invariant mass of the projectile-nucleon
*  system
*     ECMO = AIT2 + ATA2 + 2.D+00 * AMMTA * ELAB
*  Put here the correct value for Ecmo taking into account
*  Fermi motion etc.
*  Compute Ekrecl in the most favourable situation, with final
*  excitation energy zero (tvepsi) and all the energy given to the
*  projectile
*  Dex = - Tveuz, Ekres = sqrt(amnres**2 + pfrmi**2) - amnres,
*  Ekrecl = Ekres + Dex
*  **** Now the checks are done using atomic masses directly **** *
      EKREC0 = EKRECL
      TVEUZ0 = TVEUZ
      PFRMSQ = PXFRM**2 + PYFRM**2 + PZFRM**2
      EKRES  = SQRT ( ( AMMRES + 2.D+00 * TVEPSI )**2 + PFRMSQ ) -
     &         AMMRES - 2.D+00 * TVEPSI
*  -TVEUZ to reduce the excitation energy not the recoil one
* (it would have been more exact to set
*  Efrm = Efrm + Tveuz, but it is better to conserve Efrm
*  coherent with the original value)
      EKRECL = EKRES - TVEUZ + 2.D+00 * TVEPSI
      EZERO  = EFRM - V0WELL (ITJ) - EKRECL - EBNDNG (ITJ)
      ECHCK  = ELAB + EZERO
      PLABSQ = PLAB * PLAB
      PFRSCA = PXFRM * CX + PYFRM * CY + PZFRM * CZ
      ECMO   = ECHCK**2 - PLABSQ - PFRMSQ - 2.D+00 * PLAB * PFRSCA
      UMIN2  = ( AMMIT + AMMTA )**2
*  +----------------------------------------------------------------*
*  |   Check if the Hadrin call is energetically possible, else
*  |   give the total energy to cascade and excitation:
      IF ( ECMO .LT. UMIN2 .AND. NXPI .LE. 0 ) THEN
         V0WELL (ITJ) = V0WELL (ITJ) - V0EXTR
         ANOW  = ANOW  + IBAR (IT)
         KTARN = KTARN - IBAR (IT) + ICH (IT)
         ZNOW  = ZNOW  + ICH  (IT)
         KTARP = KTARP - ICH  (IT)
         TPK  = 0.D+00
         TNK  = 0.D+00
         TV   = 0.D+00
         TPNVK  = 0.D+00
         TOEFF  = ELAB - IBAR (IT) * AMMIT
         LCORIN = .TRUE.
         CALL CORSTP ( TOEFF )
         CALL GRANOR(RANGS1, RANGS2)
         TNKTE = MAX ( 0.003D+00, EKMNNU (2)
     &         * ( 0.5D+00 - 0.25D+00 * RANGS1 ) )
         TPKTE = MAX ( 0.003D+00, EKMNNU (1)
     &         * ( 0.5D+00 - 0.25D+00 * RANGS2 ) )
         GO TO 405
*  |
*  |-->-->-->-->--> go to the cascade simulation
      END IF
*  |
*  +-------------------------------------------------------------------*
      EKRECL = EKREC0
*  ******  The check for peripheral collisions is now moved here: *****
*   The old comment was:
*   " I.E.,H-COLLISION WITH N IN THE MOST OUTSIDE NUCLEUS SHELL
*     WITH THE THICKNESS OF ONE NUCLEON RADIUS
*     RATIO OF THE SHELL VOLUME TO NUCLEUS VOLUME IS ELANU3
*     THE H-N-COLL. WILL BE SAMPLED FIRSTLY WITH THE TOTAL PRIMARY
*     ENERGY ELAB, THE REMAINING ENERGY IS USED IN THE ABOVE SAMPLED
*     RATIO FOR CASCADE AND EXCITATION SWITCH VARIABLE IS IPERCO=1 "
*   The old coding was: ....
*     ELRAN  = RNDM (V)
*     ELANU3 = ANUC**(-0.6666666666666667D+00)
*     IF ( ELRAN .LE. ELANU3 .AND. ABS (ELPP-T3) .GT. 0.0001D0)
*    ..... but
*    in my opinion we must perform peripheral collisions checking
*    the ratio between the total area (r**2) and the circular
*    corona corresponding to 1 nucleon radius, so checking:
*    (now to switch off peripheral collisions we need Elanu3 = 1)
*    Furthermore no cascade and excitation energy is computed since
*    it is extremely likely that the remaining energy (not more than
*    Tveuz) will go thouroly in excitation one: anyway it is divided
*    among the 3 contributions in the usual way
*     ELRAN  = RNDM (V)
*     ELANU3 = ( ( ANUC - 1.D+00 ) / ANUC )**0.6666666666666667D+00
*    We changed again to the original idea, which can be geometrically
*    explained as a corona of dr = r0 , corresponding roughly to
*    1 nucleon. Anyway we need a much more detailed insight of this
*    problem with possibly checks with experimental data
      CALL GRNDM(RNDM,1)
      ELRAN  = RNDM (1)
      ELANU3 = 1.D+00 - 1.D+00 / ANUC**0.6666666666666667D+00
*  +-------------------------------------------------------------------*
*  |  Check for peripheral collisions: note that in this way peripheral
*  |  collisions are switched off for stopping antibaryons!!
      IF ( ELRAN .GE. ELANU3 .AND. NXPI .EQ. 0 ) THEN
         IPERCO = 1
         TOEFF  = ELAB - AMMIT
         LCORIN = .TRUE.
         CALL CORSTP ( TOEFF )
         CALL GRANOR(RANGS1, RANGS2)
         TNKTE = MAX ( 0.003D+00, EKMNNU (2)
     &         * ( 0.5D+00 - 0.25D+00 * RANGS1 ) )
         TPKTE = MAX ( 0.003D+00, EKMNNU (1)
     &         * ( 0.5D+00 - 0.25D+00 * RANGS2 ) )
         EKRECL = EKREC0
         ELPP   = ELAB
         PLPP   = PLAB
         DEX    = 0.D+00
         PSPESQ = PLABSQ + PFRMSQ + 2.D+00 * PLAB * PFRSCA
         IUMO   = 0
*  |  +----------------------------------------------------------------*
*  |  |
1500     CONTINUE
            IUMO   = IUMO   + 1
            EEX    = TVEUZ  + DEX
*  |  |  Now for iperco we are working with atomic masses
*           AMSTAR = AMNRES + EEX
            AMSTAR = AMMRES + EEX
            EKRES  = SQRT ( AMSTAR**2 + PFRMSQ ) - AMSTAR
*  |  |   Note +DEX to reduce actually the excitation energy, not
*  |  |   the recoil one (it would have been more exact to set
*  |  |   Efrm = Efrm - Dex, but it is better to conserve Efrm
*  |  |   coherent with the original value)
            EKRECL = EKRES + DEX
            EZERO  = EFRM  - V0WELL (ITJ) - EKRECL - EBNDNG (ITJ)
            ESPENT = ELPP  + EZERO
*           ELEFT  = ETTOT - ESPENT
            ELEFT  = ETTOT - ESPENT + DELEFT
            UMO2   = ESPENT**2 - PSPESQ
            DELTU2 = UMIN2 - UMO2
*  |  |  +-------------------------------------------------------------*
*  |  |  |  We need more invariant mass for the collision!!!
            IF ( DELTU2 .GT. 0.D+00 ) THEN
*  |  |  |  +----------------------------------------------------------*
*  |  |  |  |
               IF ( IUMO .GT. 3 ) THEN
                  V0WELL (ITJ) = V0WELL (ITJ) - V0EXTR
                  ANOW  = ANOW  + IBAR (IT)
                  KTARN = KTARN - IBAR (IT) + ICH (IT)
                  ZNOW  = ZNOW  + ICH  (IT)
                  KTARP = KTARP - ICH  (IT)
                  GO TO 405
               END IF
*  |  |  |  |
*  |  |  |  +----------------------------------------------------------*
               DELTEX = 0.5D+00 * DELTU2 * ELEFT / ( ESPENT * AMSTAR )
               DELTEX = MIN ( DELTEX, EEX )
               DEX    = DEX - DELTEX
            GO TO 1500
*  |  |  |
*  |  +-<|--<--<--<--<--< go to compute the new invariant mass!!
*  |  |  |
*  |  |  +-------------------------------------------------------------*
*  |  |  |
            ELSE
               TVEUZ = EEX
            END IF
*  |  |  |
*  |  |  +-------------------------------------------------------------*
*  |  |
*  |  +----------------------------------------------------------------*
*  |  we are ready for the interaction !!!!!
         IHACA = 0
*  |  +----------------------------------------------------------------*
*  |  |
1600     CONTINUE
            ITFRHD = IT
            IHACA  = IHACA + 1
            IF ( PLPP .LT. 1.D-04 ) THEN
               WRITE (LUNERR,*)
     &              ' Nucriv: iperco=1,plpp,elpp,it',
     &                PLPP,ELPP,IT
               WRITE (LUNERR,*)'      IHACA,LVMASS,ELAB,PLAB',
     &                                IHACA,LVMASS,ELAB,PLAB
            END IF
            CALL FERHAV ( IT, ELPP, PLPP, CX, CY, CZ )
*  |  |  +-------------------------------------------------------------*
*  |  |  | Check if we have reached the maximum numbers of calls
            IF ( IHACA .GT. 15 ) THEN
               GO TO 1700
*  |  |  |-->-->-->-->--> give up !!!
            END IF
*  |  |  |
*  |  |  +-------------------------------------------------------------*
*  |  | Ir is the number of secondaries produced in the
*  |  | Ferhad/Hadrin call
*  |  |  +-------------------------------------------------------------*
*  |  |  |
            IF ( IR .EQ. 0 ) THEN
               V0WELL (ITJ) = V0WELL (ITJ) + V0EXTR
               GO TO 1600
*  |  |  |
*  |  +-<|-<--<--<--<--< go to resampling the interaction
            END IF
*  |  |  |
*  |  |  +-------------------------------------------------------------*
*  |  | Of course ir=1 means that the only outgoing particle is
*  |  | still the projectile !! If we have an incoming antibaryon
*  |  | we need at least 2 particles in the final state (why???)
*  |  |  +-------------------------------------------------------------*
*  |  |  |
            IF ( IR .LT. 2 .AND. INXPI .LT. 0 ) THEN
               V0WELL (ITJ) = V0WELL (ITJ) + V0EXTR
               GO TO 1600
*  |  |  |
*  |  +-<|-<--<--<--<--< go to resampling the interaction
            END IF
*  |  |  |
*  |  |  +-------------------------------------------------------------*
1700     CONTINUE
         TVEUZ = TVEUZ + EKRECL
         TV  = TVEUZ
         TPK = 0.D+00
         TNK = 0.D+00
         GO TO 2800
*  |  |                           End loop for Ferhad/hadrin call
*  |  +----------------------------------------------------------------*
*  |    Yes peripheral collison
*  +-------------------------------------------------------------------*
*  |    No peripheral collision
      ELSE
         IPERCO = -1
      END IF
*  |
*  +-------------------------------------------------------------------*
      I1651 = 0
*  +-------------------------------------------------------------------*
*  |  It is possible to come here for resampling !!!!
1651  CONTINUE
         TVEUZ = TVEUZ0
         IRN = NP0
         DDT = CDDT * TLAB
         EMD = 0.D+00
C
C--------------------------------------------------
C*** EFFECTIVE COLLISION ENERGY CALCULATION FOR USE IN CASCADE
C*** AND EXCITATION ENERGY CALCULATION FOR ANNIHILATION
C--------------------------------------------------
C
C--------------------------------------------------
C*** I653 LOOP COUNTER INTEGER FOR BAD ENERGY CORRECTION IN ANNIHILATION
C--------------------------------------------------
         I653 = 0
*  |  +----------------------------------------------------------------*
*  |  |  It is possible to come here for resampling !!!!
 653     CONTINUE
C
C--------------------------------------------------
C*** TOEFF = EFFECTIVE KINETIC COLLISION ENERGY
C--------------------------------------------------
            I653 = I653 + 1
            TOEFF= TLAB
C
C--------------------------------------------------
C*** TOEFF ONLY FOR ANNIHILATION NE TLAB POSSIBLE, ELSE JUST SELECTION
C*** OF ENERGY FRACTIONS
C--------------------------------------------------
*  |  |  +-------------------------------------------------------------*
*  |  |  |   If we are here with a pi- or a K- this means Nxpi = 0,
*  |  |  |   ==> no stopping (annihilating) particle but Hadrin call
*  |  |  |   for pi-, K-, while for ap, an, alambdas may be that Nxpi=1
*  |  |  |   ==> no Hadrin call (more exactly, Hadrin can be called but
*  |  |  |   the particle is stopping ( annihilating ) )
            IF ( INXPI .LE. 0 ) THEN
*  |  |  |   Note, Ibar (itta) = 1 of course, Ibar (it) maybe 0 or -1
*  |  |  |   so Ibl = 0 for pi-, K- and Ibl = 1 for ap, an, alambdas
               IBL = IBAR (ITTA) * IBAR (IT)
               IBL = ABS  (IBL)
*  |  |  |   Icl = -1 for pi-, K-, ap, Icl = 0 for an and alambdas
               ICL = ICH (IT)
*  |  |  |   ( Ibar (itta) - Ibar (it) ) * Ibl = 2 for ap, an, alambdas
*  |  |  |   ( Ibar (itta) - Ibar (it) ) * Ibl = 0 for pi-, K-
*  |  |  |   ( 1 - Ibl ) * |Icl| = 0 for ap, an, alambdas
*  |  |  |   ( 1 - Ibl ) * |Icl| = 1 for K-, pi-
               EMD = AM (IT) * ( ( ( IBAR (ITTA) - IBAR(IT) ) * IBL
     &             + ( 1 - IBL ) * ABS (ICL) ) )
*  |  |  |   Finally Emd = 2 Am(it) for antibaryons and Am (it) for
*  |  |  |   pi- and K-. Other particles cannot enter the if, however
*  |  |  |   for baryons it is 0 regardless of the charge and for
*  |  |  |   positive charged mesons is Am (it): the conclusion is
*  |  |  |   that who wrote the code is a ......... . Even not
*  |  |  |   taking into account the "if" all the stuff could be:
*              EMD = AM (IT) * ( 1 - IBAR (IT) )
C
C--------------------------------------------------
C*** EMD=2*AM(IT) FOR ANNIHILATION, =AM(IT) FOR MESONS, =0ELSE
C--------------------------------------------------
               THRES = THRESR (IT)
*  |  |  | Ecms is the total available energy in the Cms frame
*  |  |  | and of course is the invariant mass of the projectile-
*  |  |  | nucleon system
               ECMS  = SQRT  (ECMO)
               IF (ECMS .LT. THRES) EMD = 0.D0
*  |  |  | This card means Emd = 0 always for pi-, K- !! and also
*  |  |  | Emd = 0 for ap, an and alambdas if Hadrin will be called
*  |  |  | (they are not stopping!!)
               IF (NXPI .LT. 1) EMD = 0.D0
C
C--------------------------------------------------
C*** EMD=0 FOR MESONS AND BARYONS,IF NO STOPPING OF PARTICLES IF H-N-CMS
C*** -ENERGY LT TABULATED THRESHOLD
C--------------------------------------------------
*  |  |  | This card means Emd = 0 always for pi-, K- (and baryons, but
*  |  |  | for them we don't enter the if)
               IF ( IBAR(IT) .GE. 0 ) EMD = 0.D0
C--------------------------------------------------
C*** ECI = TOTAL EFFECTIVE ENERGY IN H-N-CMS-SYSTEM
C--------------------------------------------------
*  |  |  | Eci is the total energy of the projectile in the frame
*  |  |  | were the target is at rest, if the "effective" invariant
*  |  |  | mass id Ecms + Emd (see kinematics inside Ferevv and
*  |  |  | Difevv )
               ECI   = 0.5D+00 * ( ( ECMS + EMD )**2 - AIT2 - ATA2 ) /
     &                 AMMTA
               TOEFF = ECI - AMMIT
*  |  |  |  +----------------------------------------------------------*
*  |  |  |  |
               IF ( TOEFF .LT. TLAB ) THEN
                  TOEFF = TLAB
               END IF
*  |  |  |  |
*  |  |  |  +----------------------------------------------------------*
            END IF
*  |  |  |
*  |  |  +-------------------------------------------------------------*
            LCORIN = .TRUE.
*  |  | Now the entry Corrnc in Corrin set the correct values for the
*  |  | excitation energy Tv and the cascade energies Tpk and Tnk
            CALL CORRNC ( TOEFF )
            CALL GRANOR(RANGS1, RANGS2)
            TNKTE = MAX ( 0.003D+00, EKMNNU (2)
     &            * ( 0.5D+00 - 0.25D+00 * RANGS1 ) )
            TPKTE = MAX ( 0.003D+00, EKMNNU (1)
     &            * ( 0.5D+00 - 0.25D+00 * RANGS2 ) )
            TV  = TVGRE0
            TPK = EINCP
            TNK = EINCN
C--------------------------------------------------
C*** ELPP=TOTAL PROJECTILE ENERGY,PLPP=ABSOLUT MOMENTUM OF THE PROJECTIL
C*** IN THE LABSYSTEM
C
C--------------------------------------------------
            T3   = AM (IT)
            TPNVK = TV + TNK + TPK
            ELPP  = ELAB - TPNVK
*  |  |  +-------------------------------------------------------------*
*  |  |  | Check for stopping pi- and K-: in my mind this point can
*  |  |  | never be reached by stopping pi-, K- since they should be
*  |  |  | already gone to 405. Anyway...
            IF ( INXPI .EQ. 0 .AND. NXPI .NE. 0 ) THEN
*  |  |  |  +----------------------------------------------------------*
*  |  |  |  |    Again in my mind Inxpi = 0 <==> Ibar (itttt) = 0,
*  |  |  |  |    anyway ....
               IF ( IBAR (ITTTT) .GE. 0 ) THEN
                  TVEUZ = 0.D+00
                  V0WELL (ITJ) = V0WELL (ITJ) - V0EXTR
                  ANOW  = ANOW  + IBAR (IT)
                  KTARN = KTARN - IBAR (IT) + ICH (IT)
                  ZNOW  = ZNOW  + ICH  (IT)
                  KTARP = KTARP - ICH  (IT)
                  GO TO 405
*  |  |  |  |-->-->-->-->--> go to the cascade simulation !!
               END IF
*  |  |  |  |
*  |  |  |  +----------------------------------------------------------*
            END IF
*  |  |  |
*  |  |  +-------------------------------------------------------------*
C
C--------------------------------------------------
C*** IF NO PRIMARY PI-,K- OR IF HADRIN FOR PI-,K- WILL BE CALLED
C*** (NO STOPPING PARTICLE CASE) GO TO 403
C--------------------------------------------------
C
C--------------------------------------------------
C*** CASE OF ENERGY CORRECTION IF LOOP IN ANNIHILATION
C--------------------------------------------------
*  |  |  +-------------------------------------------------------------*
*  |  |  |
            IF ( IBAR (IT) .GE. 0 .AND. NXPI .EQ. 0 ) THEN
               ELIM = T3 + TVEPSI
*  |  |  |
*  |  |  +-------------------------------------------------------------*
*  |  |  |
            ELSE
               ELIM = AM (13)
            END IF
*  |  |  |
*  |  |  +-------------------------------------------------------------*
 6530       CONTINUE
*  |  |  +-------------------------------------------------------------*
*  |  |  |   If Elpp < 0 for less then 10 iterations (I653), start
*  |  |  |   again with effective collision energy calculation, else
*  |  |  |   correct TPK, TNK and TV. Now the maximum number of itera-
*  |  |  |   tion has been decreased to 4, and we check Elpp < T3.
*  |  |  |   We have no enough energy or we have looped too many
*  |  |  |   times, so give up and force the "effective" energy to
*  |  |  |   be 0, reducing cascade and excitation energy
            IF ( I653 .GE. 4 .AND. ELPP .LT. ELIM ) THEN
*  |  |  |   The following cards simply take away the negative
*  |  |  |   "Elpp" to cascade and excitation in the same proportion
*  |  |  |   they were
*  |  |  |   Modified to take into account the Fermi energy:
               A65 = 1.D+00 + ( ELPP - ELIM ) / TPNVK
               TPK = TPK * A65
               TNK = TNK * A65
               TV  = TV  * A65
               TPNVK  = TPNVK * A65
               TVGRE0 = TV
               ELPP   = ELIM
*  |  |  |
*  |  |  +-------------------------------------------------------------*
*  |  |  |
            ELSE IF ( ELPP .LT. 0.D+00 ) THEN
               GO TO 653
*  |  |  |
*  |  +-<|--<--<--<--< go to resampling Tv, Tnk, Tpk
            END IF
*  |  |  |
*  |  |  +-------------------------------------------------------------*
 
*  |  |  +-------------------------------------------------------------*
*  |  |  |              If the available energy is larger than the
*  |  |  |              mass energy of the projectile check the inva-
*  |  |  |              riant mass of the system
            IF ( ELPP .GT. T3 .AND. NXPI .EQ. 0 ) THEN
               IUMO = 0
               DEX  = 0.D+00
*  |  |  |  +----------------------------------------------------------*
*  |  |  |  |
2500           CONTINUE
*  |  |  |  |   Now also here work with atomic masses
                  IUMO = IUMO + 1
                  PLPPSQ = ( ELPP + AMMIT ) * ( ELPP - AMMIT )
                  PLPP   = SQRT ( PLPPSQ )
                  PPLUS  = PLAB - PLPP
                  EPLUS  = ELAB - ELPP
                  EEX    = TVEUZ  + DEX
                  AMSTAR = AMMRES + EEX
                  EKRES  = SQRT ( AMSTAR**2 + PPLUS**2 + PFRMSQ - 2.D+00
     &                   * PFRSCA * PPLUS ) - AMSTAR
*  |  |  |  |   Note +DEX to reduce actually the excitation energy, not
*  |  |  |  |   the recoil one (it would have been more exact to set
*  |  |  |  |   Efrm = Efrm - Dex, but it is better to conserve Efrm
*  |  |  |  |   coherent with the original value)
                  EKRECL = EKRES - EPLUS + DEX
*  |  |  |  |  +-------------------------------------------------------*
*  |  |  |  |  |  Check the recoil energy:
                  IF ( EKRECL .LT. 0.D+00 .AND. IUMO .LE. 2 ) THEN
                     EKRECL = EKREC0
                     LDELTX = .FALSE.
*  |  |  |  |  |
*  |  |  |  |  +-------------------------------------------------------*
*  |  |  |  |  |
                  ELSE
                     LDELTX = .TRUE.
                  END IF
*  |  |  |  |  |
*  |  |  |  |  +-------------------------------------------------------*
                  EZERO  = EFRM  - V0WELL (ITJ) - EKRECL - EBNDNG (ITJ)
                  ESPENT = ELPP  + EZERO
                  ELEFT  = ETTOT - ESPENT
                  EELEFT = ELEFT + DELEFT
                  UMO2   = AIT2 + EZERO**2 - PFRMSQ + 2.D+00 * ( EZERO
     &                   * ELPP - PLPP * PFRSCA )
                  DELTU2 = UMIN2 - UMO2
*  |  |  |  |  +-------------------------------------------------------*
*  |  |  |  |  |  We need more invariant mass for the collision!!!
                  IF ( DELTU2 .GT. 0.D+00 ) THEN
*  |  |  |  |  |  +----------------------------------------------------*
*  |  |  |  |  |  |
                     IF ( IUMO .GT. 4 ) THEN
                        GO TO 653
*  |  |  |  |  |  |
*  |  +-<|-<|-<|-<|--<--<--<
                     END IF
*  |  |  |  |  |  |
*  |  |  |  |  |  +----------------------------------------------------*
 
*  |  |  |  |  |  +----------------------------------------------------*
*  |  |  |  |  |  |
                     IF ( LDELTX .AND. TPNVK .GT. 0.D+00 ) THEN
                        AHELP  = ESPENT / ELEFT
                        FFRAC  = ( EEX - TVEPSI ) / TPNVK
                        BHELP  = - ELPP - ELPP / PLPP * ( PFRSCA + AHELP
     &                         * ( PFRSCA - PPLUS ) )
                        DELTAE = 0.5D+00 * DELTU2 / ( BHELP + FFRAC *
     &                           AMSTAR * AHELP )
                        TMPDEL = 1.01D+00 * DELTAE
                        DELTAE = MIN ( TMPDEL, TPNVK )
                        DELTEX = FFRAC * DELTAE
                        TKOR = 1.D+00 - DELTAE / TPNVK
                        TPK  = TPK * TKOR
                        TNK  = TNK * TKOR
                        TV   = TV  * TKOR
                        TVGRE0 = TV
                        TPNVK  = TPNVK - DELTAE
*  |  |  |  |  |  |
*  |  |  |  |  |  +----------------------------------------------------*
*  |  |  |  |  |  |
                     ELSE IF ( TPNVK .GT. 0.D+00 ) THEN
                        DELTAE = 0.5D+00 * DELTU2 / ( EZERO - PFRSCA *
     &                           ELPP / PLPP )
                        TMPDEL = 1.01D+00 * DELTAE
                        DELTAE = MIN ( TMPDEL, TPNVK )
                        DELTEX = 0.D+00
                        TKOR = 1.D+00 - DELTAE / TPNVK
                        TPK  = TPK * TKOR
                        TNK  = TNK * TKOR
                        TV   = TV  * TKOR
                        TVGRE0 = TV
                        TPNVK  = TPNVK - DELTAE
*  |  |  |  |  |  |
*  |  |  |  |  |  +----------------------------------------------------*
*  |  |  |  |  |  |
                     ELSE IF ( LDELTX ) THEN
                        DELTAE = 0.D+00
                        DELTEX = 0.5D+00 * DELTU2 * ELEFT / ( ESPENT *
     &                           AMSTAR )
                        TMPDEL = 1.01D+00 * DELTEX
                        DELTEX = MIN ( TMPDEL, EEX - TVEPSI )
*  |  |  |  |  |  |
*  |  |  |  |  |  +----------------------------------------------------*
*  |  |  |  |  |  |
                     ELSE
                        DELTAE = 0.D+00
                        DELTEX = 0.D+00
                     END IF
*  |  |  |  |  |  |
*  |  |  |  |  |  +----------------------------------------------------*
*                    ELPP   = ELPP  + DELTAE, AMMIT
                     TMPAMM = AMMIT + 0.01D+00
                     ELPP   = MAX ( ELPP  + DELTAE, TMPAMM )
                     DEX    = DEX   - DELTEX
                     GO TO 2500
*  |  |  |  |  |
*  |  |  |  +-<|--<--<--<--<--< go to compute the new invariant mass!!
*  |  |  |  |  |
*  |  |  |  |  +-------------------------------------------------------*
*  |  |  |  |  |
                  ELSE
                     TVEUZ = EEX
                     GO TO 510
*  |  |  |  |  |
*  |  |  |  |  |-->-->-->-->--> invariant mass larger than the projec-
*  |  |  |  |  |                tile mass plus target nucleon mass, go
*  |  |  |  |  |                to the event simulation !!
                  END IF
*  |  |  |  |  |
*  |  |  |  |  +-------------------------------------------------------*
*  |  |  |  |
*  |  |  |  +----------------------------------------------------------*
            END IF
*  |  |  |
*  |  |  +-------------------------------------------------------------*
            IF ( IBAR (IT) .GE. 0 ) GO TO 653
*  |  |
*  |  +--<--<--<--<--< go to resampling Tv, Tnk, Tpk, Elpp < T3 and no
*  |  |                antinucleon projectile
            PLPPSQ = ( ELPP - AMMIT ) * ( ELPP + AMMIT )
*  |  |  +-------------------------------------------------------------*
*  |  |  |
            IF ( PLPPSQ .GT. 0.D+00 ) THEN
               PLPP = SQRT ( PLPPSQ )
*  |  |  |
*  |  |  +-------------------------------------------------------------*
*  |  |  |
            ELSE
               PLPP = 0.D+00
            END IF
*  |  |  |
*  |  |  +-------------------------------------------------------------*
            PPLUS  = PLAB - PLPP
            EPLUS  = ELAB - ELPP
            EEX    = TVEUZ
            AMSTAR = AMNRES + EEX
            EKRES  = SQRT ( AMSTAR**2 + PPLUS**2 + PFRMSQ - 2.D+00
     &             * PFRSCA * PPLUS ) - AMSTAR
            EKRECL = EKRES - EPLUS
*  |  |  +-------------------------------------------------------------*
*  |  |  |  Check the recoil energy:
            IF ( EKRECL .LT. EKREC0 ) THEN
               EKRECL = EKREC0
            END IF
*  |  |  |
*  |  |  +-------------------------------------------------------------*
            EZERO  = EFRM  - V0WELL (ITJ) - EKRECL - EBNDNG (ITJ)
            ESPENT = ELPP  + EZERO
            ELEFT  = ETTOT - ESPENT
            UMO2   = AIT2 + EZERO**2 - PFRMSQ + 2.D+00 * ( EZERO
     &             * ELPP - PLPP * PFRSCA )
            DELTU2 = UMIN2 - UMO2
            IF ( DELTU2 .LT. 0.D+00 ) GO TO 510
*  |  |
*  |  |-->-->-->-->--> invariant mass larger than the projec-
*  |  |                tile mass plus target nucleon mass, go
*  |  |                to the event simulation !!
*  |  |  If we are here we are dealing with an antibaryon and the total
*  |  |  available energy in the centre of mass frame is less than the
*  |  |  projectile mass plus the target nucelon mass
C--------------------------------------------------
C*** IF ELPP GREATERTHEN THE MASS OF IT, GO TO MOMENTUM LIMITATION FOR H
C*** ELSE IF NO EFFECTIVE  ENERGY WAS CALCULATED DO THE SAME
C*** ELSE IF IT IS NO ANTINUCLEON AND ELPP LT MASS OF IT, START AGAIN WI
C***   EFFECTIVE COLLISION ENERGY CALCULATION
C*** ELSE NEW ENERGY DEFINITION FOR IT AND NEW MASS DEFINITION FOR BARYO
C*** AND ANTIBARYONS
C--------------------------------------------------
*  |  |  +-------------------------------------------------------------*
*  |  |  |
            IF ( EMD .LE. 0.D+00 ) THEN
               ELPP = T3
               PLPP = 0.D+00
               GO TO 1510
*  |  |  |
*  |  |  |-->-->-->-->--> Skip mass redefinition if Emd is negative or
*  |  |  | zero. In my opinion now Emd is always = 0 for non stopping
*  |  |  | antibaryons, because of the if inside the Emd stuff and I am
*  |  |  | not sure we have really to arrive here under the minimum
*  |  |  | invariant mass...
            END IF
*  |  |  |
*  |  |  +-------------------------------------------------------------*
C
C--------------------------------------------------
C*** CASE OF ANNIHILATION,PSEUDO MASS DEFINITION
C--------------------------------------------------
            PLPP = 0.D+00
            UMO2 = ( ELPP + EZERO )**2 - PFRMSQ
            T3   = 0.5D+00 * SQRT ( UMO2 )
            IUMO = 0
            INEG = 0
*  |  |  +-------------------------------------------------------------*
*  |  |  |
5050        CONTINUE
               T3   = 0.9D+00 * T3
               T3SQ = T3 * T3
               PLPPSQ = ELPP**2 - T3SQ
*  |  |  |  +----------------------------------------------------------*
*  |  |  |  |
               IF ( PLPPSQ .GT. 0.D+00 ) THEN
                  IUMO = IUMO + 1
                  PLPP = SQRT ( PLPPSQ )
*  |  |  |  |
*  |  |  |  +----------------------------------------------------------*
*  |  |  |  |
               ELSE
                  INEG = INEG + 1
*  |  |  |  |  +-------------------------------------------------------*
*  |  |  |  |  |
                  IF ( INEG .GT. 100 ) THEN
                     WRITE (LUNOUT,*)
     &               ' Nucriv: impossible to get the pseudo-mass!!',
     &               ' Plpp always negative, Plpp,Elpp', PLPP, ELPP
                     WRITE (LUNERR,*)
     &               ' Nucriv: impossible to get the pseudo-mass!!',
     &               ' Plpp always negative, Plpp,Elpp', PLPP, ELPP
                     V0WELL (ITJ) = V0WELL (ITJ) - V0EXTR
                     ANOW  = ANOW  + IBAR (IT)
                     KTARN = KTARN - IBAR (IT) + ICH (IT)
                     ZNOW  = ZNOW  + ICH  (IT)
                     KTARP = KTARP - ICH  (IT)
                     GO TO 405
                  END IF
*  |  |  |  |  |
*  |  |  |  |  +-------------------------------------------------------*
                  GO TO 5050
               END IF
*  |  |  |  |
*  |  |  |  +----------------------------------------------------------*
               PPLUS  = PLAB - PLPP
               EPLUS  = ELAB - ELPP
               EEX    = TVEUZ
               AMSTAR = AMNRES + EEX
               EKRES  = SQRT ( AMSTAR**2 + PPLUS**2 + PFRMSQ - 2.D+00
     &                * PFRSCA * PPLUS ) - AMSTAR
               EKRECL = EKRES - EPLUS
*  |  |  |  +----------------------------------------------------------*
*  |  |  |  |  Check the recoil energy:
               IF ( EKRECL .LT. EKREC0 ) THEN
                  EKRECL = EKREC0
               END IF
*  |  |  |  |
*  |  |  |  +----------------------------------------------------------*
               EZERO  = EFRM - V0WELL (ITJ) - EKRECL - EBNDNG (ITJ)
               UMO2   = T3SQ + EZERO**2 - PFRMSQ + 2.D+00 * ( EZERO
     &                * ELPP - PLPP * PFRSCA )
               DELTU2 = 2.D+00 * T3 - UMO2
*  |  |  |  +----------------------------------------------------------*
*  |  |  |  |
               IF ( IUMO .GT. 10 ) THEN
                  WRITE (LUNERR,*)
     &            ' Nucriv: impossible to get the pseudo-mass!!',
     &              IUMO,DELTU2
                  V0WELL (ITJ) = V0WELL (ITJ) - V0EXTR
                  ANOW  = ANOW  + IBAR (IT)
                  KTARN = KTARN - IBAR (IT) + ICH (IT)
                  ZNOW  = ZNOW  + ICH  (IT)
                  KTARP = KTARP - ICH  (IT)
                  GO TO 405
               END IF
*  |  |  |  |
*  |  |  |  +----------------------------------------------------------*
               IF ( DELTU2 .GT. 0.D+00 ) GO TO 5050
*  |  |  |
*  |  |  |--<--<--<--<--< invariant mass too small !!
*  |  |  +-------------------------------------------------------------*
*  |  |  +-------------------------------------------------------------*
*  |  |  |   Check the pseudo-mass value (note that the following check
*  |  |  |   becomes useless (except for I653 > 10)
            IF ( T3 .LT. 0.3D+00 * AM (1) .AND. I653 .LT. 10 ) THEN
               I653 = MAX ( I653 + 1, 4 )
               ELIM = 0.33D+00 * AM (1) * ELPP / T3
               T3 = AM (IT)
               GO TO 6530
*  |  |  |--<--<--<--<--< pseudo-mass is too small !!
            END IF
*  |  |  |
*  |  |  +-------------------------------------------------------------*
*  |  |  +-------------------------------------------------------------*
*  |  |  | However let a warning message if it is met!!!
            IF ( T3 .LT. 0.14D+00 .AND. INXPI .LT. 0 .AND. NXPI .GT. 0 )
     &           THEN
               NXPI  = 1
               V0WELL (ITJ) = V0WELL (ITJ) - V0EXTR
               ANOW  = ANOW  + IBAR (IT)
               KTARN = KTARN - IBAR (IT) + ICH (IT)
               ZNOW  = ZNOW  + ICH  (IT)
               KTARP = KTARP - ICH  (IT)
               GO TO 405
*  |  |  |
*  |  |  |-->-->-->-->--> go to the cascade simulation
            END IF
*  |  |  |
*  |  |  +-------------------------------------------------------------*
 
*  |  |  +-------------------------------------------------------------*
*  |  |  |  Enter this if only with anti-nucleons!!!!!
            IF ( IBAR (ITTTT) .LT. 0 ) THEN
               ETDD  = 0.D0
*  |  |  |  Etkor
               ETKOR = TPNVK + AMMIT
1800           CONTINUE
*  |  |  |  +----------------------------------------------------------*
*  |  |  |  | Check if it was a stopping particle or not !!
               IF ( NXPI .LE. 0 ) THEN
                  ET = ELAB
C
C--------------------------------------------------
C*** THE PARTICLE IT (=AN) IS NOT STOPPED IN THE NUCLEUS
C--------------------------------------------------
                  ETREST = ET - ETKOR
                  IF ( ETREST .LT. 0.D0 ) ETDD = ETREST
                  TKOR = 1.D0 + ETDD / ( TPNVK + 1.D-10 )
*  |  |  |  |  +-------------------------------------------------------*
*  |  |  |  |  |   Treat it as a stopped particle
                  IF ( TKOR .LT. 1.D-6 ) THEN
                     NXPI = 1
                     GO TO 1800
*  |  |  |  |  |
*  |  |  |  |  +-------------------------------------------------------*
*  |  |  |  |  |        update tv, tnk, tpk to account for etdd < 0
                  ELSE
                     TV  = TV  * TKOR
                     TPK = TPK * TKOR
                     TNK = TNK * TKOR
                     TVGRE0 = TV
                     TPNVK = TPK + TNK + TV
                     ELPP= ET - TPNVK
                     GO TO 510
                  END IF
*  |  |  |  |  |
*  |  |  |  |  +-------------------------------------------------------*
*  |  |  |  |
*  |  |  |  +----------------------------------------------------------*
*  |  |  |  |       It is a stopping antibaryon !!
               ELSE
                  ETDD  = 0.D0
C
C--------------------------------------------------
C*** THE PARTICLE IT (=AN) IS STOPPED IN THE NUCLEUS
C--------------------------------------------------
                  ET = ELAB + AMMTA
                  ETREST = ET - ETKOR
                  IF ( ETREST .LT. 0.D0 ) ETDD = ETREST
                  TKOR = 1.D+00 + ETDD / TPNVK
                  IF ( TKOR .LT. 1.D-6 ) GO TO 653
*  |  |  |  |
*  |  +-<|-<|--<--<--< go to resampling Tv, Tnk, Tpk
                  TV  = TV  * TKOR
                  TPK = TPK * TKOR
                  TNK = TNK * TKOR
                  TVGRE0 = TV
                  TPNVK = TPK + TNK + TV
                  ELPP = ELAB - TPNVK
                  PLPP = SQRT ( ELPP**2 - T3SQ )
               END IF
*  |  |  |  |
*  |  |  |  +----------------------------------------------------------*
            END IF
*  |  |  |
*  |  |  +-------------------------------------------------------------*
*  |  | end of the "653" loop for bad energy correction in annihilation
*  |  +----------------------------------------------------------------*
C
C--------------------------------------------------
C*** NEW MASS DEFINITION OF BARYONS AND ANTIBARYONS
C*** (ONLY FOR REMAINING EFFECTIVE CM-ENERGIES,LOWER THAN THE TOTAL
C*** 1.88GEV-THRESHOLD CM-ENERGY IN ANNIHILATION)
C--------------------------------------------------
*  |  +----------------------------------------------------------------*
*  |  |  First pass through mass redefinition
         IF ( .NOT. LVMASS ) THEN
            LVMASS = .TRUE.
*  |  |  +-------------------------------------------------------------*
*  |  |  |  Loop over masses to be redefined
            DO 1511 IAM=1,15
               JAM        = IAMC (IAM)
               AMHH (IAM) = AM   (JAM)
               IF ( IBAR (JAM) .NE. 0 ) AM (JAM) = MIN  ( T3, AM (JAM) )
1511        CONTINUE
*  |  |  |
*  |  |  +-------------------------------------------------------------*
*  |  |
*  |  +----------------------------------------------------------------*
*  |  |  Masses already redefined
         ELSE
*  |  |  +-------------------------------------------------------------*
*  |  |  |  Loop over masses to be redefined
            DO 15111 IAM=1,15
               JAM        = IAMC (IAM)
               IF ( IBAR (JAM) .NE. 0 ) AM (JAM) = MIN ( T3, AMHH (IAM))
15111       CONTINUE
*  |  |  |
*  |  |  +-------------------------------------------------------------*
         END IF
*  |  |
*  |  +----------------------------------------------------------------*
*  | We will arrive here also with T3 = AM (IT) and ELPP < T3 from
*  | the if on EMD =< 0, it is not clear if it is correct!!
1510     CONTINUE
*  |  +----------------------------------------------------------------*
*  |  |                             Loop for Ferhad/hadrin call
510      CONTINUE
C
C--------------------------------------------------
C*** MOMENTUM LIMITATION FOR HADRIN IN FERMI-MOM.-VERSION
C--------------------------------------------------
            PLABOU=15.D0
*  |  |  +-------------------------------------------------------------*
*  |  |  |  This check seems to me to be useless, since the maximum
*  |  |  |  allowed momentum is 10 GeV/c .... anyway get a warning
*  |  |  |  message if this condition is met
            IF ( PLPP .GE. PLABOU ) THEN
               WRITE ( LUNERR, * )' Nucriv: momentum limitation met!!',
     &                              PLPP, PLABOU
               APL  = PLPP - PLABOU
               ELPP = SQRT (PLABOU**2 + AM(IT)**2)
               TPK  = TPK + APL*.3D0
               TNK  = TNK + APL*.3D0
               TV   = TV  + APL*.4D0
               IFORBI = IFORBI + 1
               IF ( IFORBI .GT. 10 ) GO TO 405
            END IF
*  |  |  |
*  |  |  +-------------------------------------------------------------*
            EELAB = ELPP
            PPLAB = PLPP
            IHACA = IHACA + 1
            ITFRHD = IT
            IF ( PPLAB .LT. 1.D-04 ) THEN
               WRITE (LUNERR,*)
     &              ' Nucriv: iperco=0,eelab,pplab,it',
     &                EELAB,PPLAB,IT
               WRITE (LUNERR,*)'      IHACA,LVMASS,ELAB,PLAB',
     &                                IHACA,LVMASS,ELAB,PLAB
            END IF
            CALL FERHAV ( IT, EELAB, PPLAB, CX, CY, CZ )
C
C--------------------------------------------------
C***
C***
C*** HADRON NUCLEON COLLISION SIMULATION
C*** EELAB=LABENERGY OF H,PPLAB=ABSOLUT MOMENTUM OF H, CX,CY,CZ DIRECTIO
C*** COSINES OF H,IN THE LAB SYSTEM, ITTA TARGET NUCLEON INDEX,
C*** ANUC,ZNUC NUCLEON AND PROTON NUMBERS
C***
C*** REGARDING FERMI MOMENTUM OF THE NUCLEON
C***
C--------------------------------------------------
C--------------------------------------------------
C***
C***
C*** IHACA TIMES CALLED H-N-COLL., IF IN FINAL STATE NO PARTICLE,CALL
C*** AGAIN,
C*** IF IN ANNIHILATION OF A STOPPING NUCLEON,   IN THE FINAL
C*** STATE IS A BARYON, CALL AGAIN H-N-COLL.
C*** (NOT MORE THAN 30 TIMES, ELSE TAKE THIS STATE)
C***
C***
C--------------------------------------------------
*  |  |  +-------------------------------------------------------------*
*  |  |  | Check if we have reached the maximum numbers of calls (which
*  |  |  | has been reduced to 15 .... 30 seems to be a bit too large
*  |  |  | but there were two "IHACA = IHACA + 1" ....
            IF ( IHACA .GT. 15 ) THEN
              GO TO 88
*  |  |  |-->-->-->-->--> give up !!!
            END IF
*  |  |  |
*  |  |  +-------------------------------------------------------------*
*  |  | Ir is the number of secondaries produced in the Ferhad/Hadrin
*  |  | call
*  |  |  +-------------------------------------------------------------*
*  |  |  |
            IF (IR .EQ. 0) THEN
               V0WELL (ITJ) = V0WELL (ITJ) + V0EXTR
               GO TO 510
*  |  |  |
*  +-<|-<|-<--<--<--<--< go to resampling the ferhad/hadrin interaction
            END IF
*  |  |  |
*  |  |  +-------------------------------------------------------------*
*  |  | Check that for annihilated antibaryons we have no (anti)baryon
*  |  | in the final state: if yes resample !
*  |  |  +-------------------------------------------------------------*
*  |  |  |
            IF ( INXPI .LT. 0 .AND. NXPI .GT. 0 ) THEN
               IBAROT = 0
               DO 520 IRZ = 1, IR
                  IBAROT = IBAROT + ABS (IBAR(ITR(IRZ)))
  520          CONTINUE
               IF ( IBAROT .GT. 0 ) THEN
                  V0WELL (ITJ) = V0WELL (ITJ) + V0EXTR
                  GO TO 510
               END IF
*  |  |  |
*  +-<|-<|-<--<--<--<--< go to resampling the ferhad/hadrin interaction
            END IF
*  |  |  |
*  |  |  +-------------------------------------------------------------*
C
C--------------------------------------------------
C*** IF NO ANNIHILATION INTO MESONS FOR AN+N
C--------------------------------------------------
*  |  | Of course ir=1 means that the only outgoing particle is still
*  |  | the projectile !! If we have an incoming antibaryon
*  |  | we need at least 2 particles in the final state (why????, if
*  |  | we have no annihilation for example a charge exchange recation
*  |  | can give ir=2)
*  |  |  +-------------------------------------------------------------*
*  |  |  |
            IF ( IR .LT. 2 .AND. INXPI .LT. 0 ) THEN
               V0WELL (ITJ) = V0WELL (ITJ) + V0EXTR
               GO TO 510
*  |  |  |
*  +-<|-<|-<--<--<--<--< go to resampling the ferhad/hadrin interaction
            END IF
*  |  |  |
*  |  |  +-------------------------------------------------------------*
88       CONTINUE
*  |  |                                End loop for Ferhad/hadrin call
*  |  +----------------------------------------------------------------*
 
*  |  +----------------------------------------------------------------*
*  |  |  Redefinition of masses, if they were changed for annihilation
*  |  |  beyond 1.88 GeV - threshold
         IF ( LVMASS ) THEN
*  |  |  +-------------------------------------------------------------*
*  |  |  |
            DO 1512 IAM=1,15
               AM (IAMC(IAM)) = AMHH (IAM)
1512        CONTINUE
*  |  |  |
*  |  |  +-------------------------------------------------------------*
         END IF
*  |  |
*  |  +----------------------------------------------------------------*
C
C--------------------------------------------------
C*** CORRECTION OF CASCADE AND EXCITATION ENERGIES FOR FERMI-MOMENTUM
C*** IF ENERGY CONSERVATION ALLOWS NO H-N-COLLISION, GO TO THE START
C*** POINT OF EVENT SIMUL.
C*** ELSE ENERGY CORRECTION BY VARIABLE EKIKOR
C--------------------------------------------------
C--------------------------------------------------
2800     CONTINUE
         ECHCK  = ETTOT
         PXCHCK = PXTTOT
         PYCHCK = PYTTOT
         PZCHCK = PZTTOT
*  |  +----------------------------------------------------------------*
*  |  |
         DO 3000 IPART = 1, IR
            ECHCK  = ECHCK  - EL (IPART)
            PXCHCK = PXCHCK - PL (IPART) * CXR (IPART)
            PYCHCK = PYCHCK - PL (IPART) * CYR (IPART)
            PZCHCK = PZCHCK - PL (IPART) * CZR (IPART)
3000     CONTINUE
*  |  |
*  |  +----------------------------------------------------------------*
         UMO2 = (ECHCK + DELEFT)**2 - PXCHCK**2 - PYCHCK**2 - PZCHCK**2
*  |  +----------------------------------------------------------------*
*  |  |    Now moved to atomic masses !!!
         IF ( IR .GT. 1 ) THEN
            AMMRE2 = AMMRES**2
*  |  |
*  |  +----------------------------------------------------------------*
*  |  |    It is supposed that the only emitted particle is still the
*  |  |    projectile
         ELSE
            AMMRES = AMMTAR
            AMNRES = AMNTAR
            AMMRE2 = AMMTAR**2
         END IF
*  |  |
*  |  +----------------------------------------------------------------*
*  |  +----------------------------------------------------------------*
*  |  |
         IF ( UMO2 .LT. AMMRE2 ) THEN
            I1651 = I1651 + 1
            IF ( I1651 .GT. 4 .OR. IHACA .GT. 10 ) THEN
               LRESMP = .TRUE.
               RETURN
            END IF
            V0WELL (ITJ) = V0WELL (ITJ) + V0EXTR
            GO TO 1651
*  |  |
*  +-<|-<--<--<--<--<--< go to the beginning of the event simulation!
         END IF
*  |  |
*  |  +----------------------------------------------------------------*
*  |              End of the event simulation: finally we got it
*  +-------------------------------------------------------------------*
*   Update Tv, Tnk, Tpk taking into account Ekikor: simply ekikor is
*   added (with its sign) leaving unchanged the ratios
*   Again this stuff is now obsolete, we make the same from the energy
*   balance
*  +-------------------------------------------------------------------*
*  |
      IF ( IPERCO .LT. 0 ) THEN
         ELRES  = ECHCK - AMNRES
         ELRES0 = TVEUZ + TPNVK
         AA65 = ELRES / ELRES0
*  | Tpnvk records the sum of the energy available for cascade protons,
*  | neutrons and excitation
*  |  +----------------------------------------------------------------*
*  |  |
         IF ( TPNVK .LE. 0.D+00 ) THEN
            TPK = 0.4D+00 * ELRES
            TNK = 0.5D+00 * ELRES
            TV  = 0.1D+00 * ELRES
            TVGRE0 = 0.5D+00 * TV
            TVEUZ  = 0.5D+00 * TV
*  |  |
*  |  +----------------------------------------------------------------*
*  |  |
         ELSE
            TPK = TPK * AA65
            TNK = TNK * AA65
            TVGRE0 = TVGRE0 * AA65
            TVEUZ  = TVEUZ  * AA65
            TV = TVGRE0 + TVEUZ
         END IF
*  |  |
*  |  +----------------------------------------------------------------*
*  |
*  +-------------------------------------------------------------------*
*  |                           peripheral collision, no cascade
      ELSE
         ELRES0 = TVEUZ
         ELRES  = ECHCK - AMNRES
         AA65 = ELRES / ELRES0
         TPK = 0.D+00
         TNK = 0.D+00
         TVGRE0 = 0.D+00
         TVEUZ  = TVEUZ * AA65
         TV = TVEUZ
      END IF
*  |
*  +-------------------------------------------------------------------*
C
C--------------------------------------------------
C
C*** REDEFINE THE ENERGIES OF SAMPLED PARTICLES IN
C*** CASE OF CHANGED MASSES
C
C--------------------------------------------------
*  +-------------------------------------------------------------------*
*  |
      IF ( LVMASS ) THEN
         MESON=0
*  |  +----------------------------------------------------------------*
*  |  |
         DO 2510 IRZ=1,IR
            ITRZ   = MIN (ITR(IRZ),30)
            IBAROT = JAMC(ITRZ)
*  |  |  +-------------------------------------------------------------*
*  |  |  |   Check if we need to redefine the energy
            IF ( ABS (IBAROT) .GE. 1 ) THEN
*  |  |  |   This energy redefinition can be very dangerous for any
*  |  |  |   energy balance, since we retain momentum and not energy!!
               ECHCK = ECHCK + EL (IRZ)
               EL (IRZ) = SQRT ( PL (IRZ)**2 + AM (ITRZ)**2 )
               ECHCK = ECHCK - EL (IRZ)
*  |  |  |
*  |  |  +-------------------------------------------------------------*
*  |  |  |
            ELSE
                MESON = MESON + 1
            END IF
*  |  |  |
*  |  |  +-------------------------------------------------------------*
2510     CONTINUE
*  |  |
*  |  +----------------------------------------------------------------*
      END IF
*  |
*  +-------------------------------------------------------------------*
 
*  +-------------------------------------------------------------------*
*  |                           No particle in the secondary stack
      IF ( IR .EQ. 0 ) THEN
         ANOW  = ANOW  + IBAR (IT)
         KTARN = KTARN - IBAR (IT) + ICH (IT)
         ZNOW  = ZNOW  + ICH  (IT)
         KTARP = KTARP - ICH  (IT)
         GO TO 405
*  |
*  +-->-->-->-->--> go to the cascade simulation
*  |
*  +-------------------------------------------------------------------*
*  |  Only one particle in the final state, it should be the projectile
*  |  since it is assumed that no nucleon number correction occurs
*  |  ( see below ): give random angle like for cascade particles to the
*  |  single secondary from Hadrin
      ELSE IF ( IR .EQ. 1 ) THEN
         AMMRES = AMMTAR
         AMNRES = AMNTAR
         TKI = EL(1) - AM (ITR(1))
         IF(TKI .LE. 1.D-20 ) GO TO 1221
         ADE=0.12D0*(1.D0+0.003D0*BBTAR)/TKI
         DEX=EXP(-1.57D0*1.57D0/ADE)
         AN1=(1.D0-DEX)*ADE*.5D0
         AN2=DEX*1.57D0
         AN=AN1+AN2
         AN1=AN1/AN
         CALL GRNDM(RNDM,1)
         IF (RNDM(1) .GT. AN1) GO TO 1222
 1223    CONTINUE
C
C--------------------------------------------------
C*** TETA ANGLE DETERMINATION (SINGLE PART.CASE), = DE
C--------------------------------------------------
         CALL GRNDM(RNDM,1)
         DE=SQRT(-ADE*LOG(1.D0-RNDM(1)*(1.D0-DEX)))
         IF(DE.GT.1.57D0) GO TO 1223
         GOTO 1224
 1222    CONTINUE
         CALL GRNDM(RNDM,1)
         DE=-RNDM(1)
         DE=ATAN2(SQRT(ONEONE-DE**2),DE)
 1224    CONTINUE
         CALL COSI(SFE,CFE)
C
C--------------------------------------------------
C*** COS PHI, SIN PHI DETERMINATION(S.P.CASE)
C--------------------------------------------------
         SID=SIN(DE)
         COD=COS(DE)
         CALL TTRANS (CXR(1), CYR(1), CZR(1), COD, SID, SFE, CFE,
     &                CCXR, CCYR, CCZR)
C
C--------------------------------------------------
C*** TURNING OF ANGLE S INTO THE LABSYSTEM (S.P.CASE)
C--------------------------------------------------
         CXR(1)=CCXR
         CYR(1)=CCYR
         CZR(1)=CCZR
 1221    CONTINUE
*  |
*  +-------------------------------------------------------------------+
*  |                           Two or more particles produced
      ELSE
         I1 = ITJ - 1
         A1 = I1
         ANOW = ANOW - 1.D0
         ZNOW = ZNOW - 1.D0 + A1
         KTARP = KTARP + 1 - I1
         KTARN = KTARN + I1
      END IF
*  |
*  +-------------------------------------------------------------------+
      IRN   = NP0
*  +-------------------------------------------------------------------*
*  |
      DO 10 I = 1,IR
C
C--------------------------------------------------
C*** STORAGE OF H-N-COLL. PRODUCTS INTO NUCLEUS FINAL STATE C./FINUC/
C*** CUT OFF FOR KINETIC ENERGIES LESS THEN ETHR, PUSH THEM INTO EX.EN.
C*** TV
C--------------------------------------------------
         T1 = EL  (I)
         I1 = ITR (I)
*  |  +----------------------------------------------------------------*
*  |  |
         IF ( T1 .LE. ETHR + AM (I1) ) THEN
            ANOW  = ANOW  + IBAR (I1)
            KTARN = KTARN - IBAR (I1) + ICH (I1)
            ZNOW  = ZNOW  + ICH (I1)
            KTARP = KTARP - ICH  (I1)
            AMMRES = ANOW * AMUAMU + 1.D-03 * FKENER ( ANOW, ZNOW )
            AMNRES = AMMRES - ZNOW * AMELEC + ELBNDE ( NINT (ZNOW) )
            EAVAIL = ECHCK - AMNRES + T1
            TPNVK  = TNK + TPK + TV
*  |  |  +-------------------------------------------------------------*
*  |  |  |
            IF ( TPNVK .GT. 0.D+00 ) THEN
               TKOR   = EAVAIL / TPNVK
               TPK = TPK * TKOR
               TNK = TNK * TKOR
               TV  = TV  * TKOR
               TVEUZ  = TVEUZ  * TKOR
               TVGRE0 = TVGRE0 * TKOR
*  |  |  |
*  |  |  +-------------------------------------------------------------*
*  |  |  |
            ELSE
               TPK = 0.4444444444444444D+00 * EAVAIL
               TNK = 0.5555555555555556D+00 * EAVAIL
            END IF
*  |  |  |
*  |  |  +-------------------------------------------------------------*
*  |  |
*  |  +----------------------------------------------------------------*
*  |  |                                  particle acceptable
         ELSE
            IRN = IRN+1
C
C--------------------------------------------------
C*** STORE THE PARTICLE CHARACTERISTICS INTO C./FINUC/, THEN TAKE NEXT
C*** PARTICLE
C--------------------------------------------------
            ITRN(IRN) = I1
            CXRN(IRN) = CXR(I)
            CYRN(IRN) = CYR(I)
            CZRN(IRN) = CZR(I)
            ELR(IRN)  = EL(I)
            PLR(IRN)  = PL(I)
            IBAT = IBAT+1
*  |  |   updating the conservation counters in common balanc
            ENUCR  = ENUCR  + ELR(IRN)
            PXNUCR = PXNUCR + PLR(IRN)*CXRN(IRN)
            PYNUCR = PYNUCR + PLR(IRN)*CYRN(IRN)
            PZNUCR = PZNUCR + PLR(IRN)*CZRN(IRN)
            ICNUCR = ICNUCR + ICH(ITRN(IRN))
            IBNUCR = IBNUCR + IBAR(ITRN(IRN))
            IRNTOT = IRNTOT + 1
         END IF
*  |  |
*  |  +----------------------------------------------------------------*
  10  CONTINUE
*  |
*  +-------------------------------------------------------------------*
 
*  +-------------------------------------------------------------------*
*  |
      IF ( IRN - NP0 - IGREYP - IGREYN .EQ. 1 ) THEN
*  |  The following balance is ok if and only if both the incident and
*  |  the emitted particles were nucleons (proton or neutron), or the
*  |  emitted particle is equal to the incident one
*  |  +----------------------------------------------------------------*
*  |  |
         IF ( ITTTT .NE. ITRN (NP0+1) ) THEN
            LPRONU = ITTTT .EQ. 1 .OR. ITTTT .EQ. 8
            LEMINU = ITRN (NP0+1) .EQ. 1 .OR. ITRN (NP0+1) .EQ. 8
*  |  |  +-------------------------------------------------------------*
*  |  |  |
            IF ( .NOT. LPRONU .OR. .NOT. LEMINU ) THEN
               ICNOW = NINT (ZNOW)
               IBNOW = NINT (ANOW)
               ICCC = ICHTAR + ICH  (ITTTT) - ICNUCR - ICINTR
               IBBB = IBTAR  + IBAR (ITTTT) - IBNUCR - IBINTR
*  |  |  |  +----------------------------------------------------------*
*  |  |  |  |
               IF ( IBBB .NE. IBNOW .AND. LEMINU .AND. IBAR (ITTTT)
     &              .EQ. 0 ) THEN
*  |  |  |  |  +-------------------------------------------------------*
*  |  |  |  |  |
                  IF ( ICCC .GT. ICNOW ) THEN
                     KTARP = KTARP - 1
                     KTARN = KTARN + 2
                     ANOW  = ANOW - 1.D+00
                     ZNOW  = ZNOW + 1.D+00
*  |  |  |  |  |
*  |  |  |  |  +-------------------------------------------------------*
*  |  |  |  |  |
                  ELSE IF (ICCC .LT. ICNOW) THEN
                     KTARP = KTARP + 1
                     KTARN = KTARN
                     ANOW  = ANOW - 1.D+00
                     ZNOW  = ZNOW - 1.D+00
*  |  |  |  |  |
*  |  |  |  |  +-------------------------------------------------------*
*  |  |  |  |  |
                  ELSE
                     KTARN = KTARN + 1
                     ANOW  = ANOW - 1.D+00
                  END IF
*  |  |  |  |  |
*  |  |  |  |  +-------------------------------------------------------*
*  |  |  |  |
*  |  |  |  +----------------------------------------------------------*
*  |  |  |  |
               ELSE IF ( ICCC .NE. ICNOW .AND. LEMINU .AND. IBAR (ITTTT)
     &                   .EQ. 0 ) THEN
*  |  |  |  |  +-------------------------------------------------------*
*  |  |  |  |  |
                  IF ( ICCC .GT. ICNOW ) THEN
                     KTARP = KTARP - 1
                     KTARN = KTARN + 1
                     ZNOW  = ZNOW + 1.D+00
*  |  |  |  |  |
*  |  |  |  |  +-------------------------------------------------------*
*  |  |  |  |  |
                  ELSE
                     KTARP = KTARP + 1
                     KTARN = KTARN - 1
                     ZNOW  = ZNOW - 1.D+00
                  END IF
*  |  |  |  |  |
*  |  |  |  |  +-------------------------------------------------------*
               END IF
*  |  |  |  |
*  |  |  |  +----------------------------------------------------------*
            END IF
*  |  |  |
*  |  |  +-------------------------------------------------------------*
         END IF
*  |  |
*  |  +----------------------------------------------------------------*
      END IF
*  |
*  +-------------------------------------------------------------------*
      GO TO 406
*  +-->-->-->-->--> go to the cascade simulation!!!
C
C--------------------------------------------------
C*** IF ITTTT IS A STOPPED NEGATIVE MESON,(BUT ALSO FOR ANNIHILATION)
C*** DISTRIBUTE THE AVAILABLE TOTAL ENERGY OF THE PARTICLE TO THE CASCAD
C*** AND EXCITATION
C--------------------------------------------------
*  +-------------------------------------------------------------------*
*  | Zero secondary part. case (h-n-coll.), stopping mesons or anyway
*  | all situations not allowing a call to Hadrin
405   CONTINUE
         V0WELL (ITJ) = V0OLD
         IRN = NP0
         AMMRES = AMUAMU * ANOW + 1.D-03 * FKENER ( ANOW, ZNOW )
         AMNRES = AMMRES - ZNOW * AMELEC + ELBNDE ( NINT (ZNOW) )
         TPNVK  = TPK + TNK + TV
         EAVAIL = ETTOT - AMNRES - TPNVK
*  |  +----------------------------------------------------------------*
*  |  |
         IF ( TPNVK .GT. 0.D+00 ) THEN
            TKOR = 1.D+00 + EAVAIL / TPNVK
            TPK  = TPK * TKOR
            TNK  = TNK * TKOR
            TV     = TV * TKOR
            TVGRE0 = TV
            TVEUZ  = 0.D+00
*  |  |
*  |  +----------------------------------------------------------------*
*  |  |
         ELSE
            TPK  = 0.4D+00 * EAVAIL
            TNK  = 0.5D+00 * EAVAIL
            TVGRE0 = 0.1D+00 * EAVAIL
            TV     = TVGRE0
            TVEUZ  = 0.D+00
         END IF
*  |  |
*  |  +----------------------------------------------------------------*
         ELPP = 0.0D+00
*  |
*  +-------------------------------------------------------------------*
*  Now the cascade simulation!!!
406   CONTINUE
      EINCP = 0.D+00
      EINCN = 0.D+00
      IU  = IRN - NP0
C
C--------------------------------------------------
C
C   SIMULATION OF CASCADE NEUTRONS
C
C--------------------------------------------------
C
C--------------------------------------------------
C*** NUCLEON NUMBER AND CASCADE ENERGY(VERSUS CUT OFF)-TEST
C*** GIVE THE ENERGY TO  THE EXCITATION , IF TEST DEMANDS THIS
C--------------------------------------------------
*  +-------------------------------------------------------------------*
*  |
      IF ((TPK.LE.TPKTE).OR.(ZNOW.LT.0.5D0)) THEN
         TV=TV+TPK
         TPK = 0.D+00
         LPCASC = .FALSE.
*  |
*  +-------------------------------------------------------------------*
*  |
      ELSE
         LPCASC = .TRUE.
      END IF
*  |
*  +-------------------------------------------------------------------*
*  +-------------------------------------------------------------------*
*  |
      IF ((TNK.LE.TNKTE).OR.(ANOW-ZNOW.LT.0.5D0)) THEN
         TV=TV+TNK
         TNK = 0.D+00
         LNCASC = .FALSE.
*  |
*  +-------------------------------------------------------------------*
*  |
      ELSE
         LNCASC = .TRUE.
      END IF
*  |
*  +-------------------------------------------------------------------*
      KNREJE = 0
      KPREJE = 0
*  +-------------------------------------------------------------------*
*  |    Cascade nucleons selection !!!!
1900  CONTINUE
*  |  +----------------------------------------------------------------*
*  |  |
         IF ( LPCASC ) THEN
            ZSAMP = ZNOW
            IF ( ZSAMP .LE. 0.5D+00 ) THEN
               LPCASC = .FALSE.
            END IF
*  |  |
*  |  +----------------------------------------------------------------*
*  |  |
         ELSE
            ZSAMP = 0.D+00
         END IF
*  |  |
*  |  +----------------------------------------------------------------*
*  |  +----------------------------------------------------------------*
*  |  |
         IF ( LNCASC ) THEN
            ANSMP = ANOW - ZNOW
            IF ( ANSMP .LE. 0.5D+00 ) THEN
               LNCASC = .FALSE.
            END IF
*  |  |
*  |  +----------------------------------------------------------------*
*  |  |
         ELSE
            ANSMP = 0.D+00
         END IF
*  |  |
*  |  +----------------------------------------------------------------*
         IF ( ( .NOT. LPCASC ) .AND. ( .NOT. LNCASC ) ) GO TO 2000
         IF ( NINT (ANOW) .LE. 0 ) THEN
            WRITE(LUNERR,*)' Nucriv: <<<< ANOW < 0 >>>>',
     &      ANOW,ZNOW,LPCASC,LNCASC
            LPCASC = .FALSE.
            LNCASC = .FALSE.
            GO TO 2000
         END IF
         CALL GRNDM(RNDM,1)
         IF ( RNDM (1) .LT. ANSMP / ( ANSMP + ZSAMP ) ) THEN
            GO TO 20
         ELSE
            GO TO 30
         END IF
   20    CONTINUE
         T2=AM(8)
C
C--------------------------------------------------
C*** NEUTRON NUMBER AND NEUTRON CASCADE ENERGY TEST
C--------------------------------------------------
         IF ( ANSMP .LT. 0.5D+00 .OR. TNK .LE. TNKTE ) THEN
            LNCASC = .FALSE.
            TV  = TV + TNK
            TNK = 0.D+00
            GO TO 1900
         END IF
C
C--------------------------------------------------
C*** SAMPLING OF KINETIC ENERGY TN AND TETA ANGLE DN OF THE EMITTED NEU
C--------------------------------------------------
         CALL RBKEKV ( 2, EXSOP, TOEFF, BBTAR, TKIN, TSTRCK,
     &                 PSTRCK, ANOW, TKRECL, COD, SID )
         TN  = TKIN
         TNK = TNK-TN
C
C--------------------------------------------------
C*** TEST: IF THE REMAINING ENERGY AND CORRESP.NUCLEON NUMBER ALLOWS FUR
C*** THER CORRESP. NUCLEON EMISSION, ELSE GIVE THE  THE ENERGY OR TAKE
C*** THE MISSING FROM THE EXCITATION
C*** IF FOR LAST CHOSEN NUCLEON TN>REMAINING REST, TAKE IN 90% THE ENER
C*** DIFFERENCE FROM EXCIT.ENERGY, IN 10% OR IF TV WOULD BE LESS THEN 0.
C***  ADD TN TO TV (It does the opposite!!! A. Ferrari)
C--------------------------------------------------
         IF ( ANSMP .LT. 1.5D+00 .OR. TNK .LE. ETHR  ) THEN
            LNCASC = .FALSE.
            CALL GRNDM(RNDM,1)
            VT = RNDM (1)
            THK=TV+TNK
            IF (VT .GT. 0.9D+00 .OR. THK .LT. 0.D+00) THEN
               TV  = TV + TNK + TN
               TNK = 0.D+00
               GO TO 1900
            END IF
            TNK=0.D0
            TV=THK
         END IF
C
C--------------------------------------------------
C*** CORRESP. NUCLEON NUMBER COUNTING
C--------------------------------------------------
         IF ( TSTRCK .LT. ETHR ) THEN
            TV = TV + TN
            GO TO 1900
*  |  |
*  +-<|--<--<--<--<--< Try to select another nucleon!!
         END IF
         ANOW   = ANOW - 1.D+00
         FRSURV = ANOW / BBTAR
C
C--------------------------------------------------
C*** TAKE THE EMITTED CORRESP. NUCLEON INTO THE FINAL PARTICLE TABLE
C*** C./FINUC/ AFTER CHOICE OF PHI AND TURNING INTO THE LAB SYSTEM
C--------------------------------------------------
         IRN = IRN + 1
         CALL COSI ( SFE, CFE )
*  |  +----------------------------------------------------------------*
*  |  |
         IF ( ANOW .LT. 30.D+00 .OR. FRSURV .LT. 0.66D+00 ) THEN
            PXLAST = PXTTOT - PXINTR - PXNUCR
            PYLAST = PYTTOT - PYINTR - PYNUCR
            PZLAST = PZTTOT - PZINTR - PZNUCR
            PPLAS2 = PXLAST**2 + PYLAST**2 + PZLAST**2
*  |  |  +-------------------------------------------------------------*
*  |  |  |
            IF ( PPLAS2 .GT. 1.D-6 * PTTOT**2 ) THEN
               PPLAST = SQRT ( PPLAS2 )
               CXXINC = PXLAST / PPLAST
               CYYINC = PYLAST / PPLAST
               CZZINC = PZLAST / PPLAST
               CDMIN = MIN ( ONEONE, PSTRCK / PPLAST ) * ( 1.D+00 -
     &                       0.8D+00 * FRSURV )
               COD = 1.D+00 - 0.5D+00 * ( 1.D+00 - COD ) * ( 1.D+00 -
     &               CDMIN )
               COD = MIN ( COD, ONEONE )
               SID = SQRT ( 1.D+00 - COD**2 )
*  |  |  |
*  |  |  +-------------------------------------------------------------*
*  |  |  |
            ELSE
               CXXINC = CX
               CYYINC = CY
               CZZINC = CZ
            END IF
*  |  |  |
*  |  |  +-------------------------------------------------------------*
            CALL TTRANS ( CXXINC, CYYINC, CZZINC, COD, SID, SFE, CFE,
     &                    CXRN (IRN), CYRN (IRN), CZRN (IRN) )
*  |  |  +-------------------------------------------------------------*
*  |  |  |
            IF ( NINT (ANOW) .GT. 0 ) THEN
               AMNRES = AMUAMU * ANOW - ZNOW * AMELEC + 1.D-03 * FKENER
     &                ( ANOW, ZNOW ) + ELBNDE ( NINT ( ZNOW ) )
*  |  |  |
*  |  |  +-------------------------------------------------------------*
*  |  |  |
            ELSE
               AMNRES = 0.D+00
            END IF
*  |  |  |
*  |  |  +-------------------------------------------------------------*
            AMNRE2 = AMNRES  * AMNRES
            ELEFT  = ETTOT - ENUCR - EINTR - TSTRCK - T2
*  |  | Update TV !!!
            TV     = ELEFT - TPK - TNK - AMNRES
            IUMO   = 0
*  |  |  +-------------------------------------------------------------*
*  |  |  |
2020        CONTINUE
               PXLEFT = PXLAST - PSTRCK * CXRN (IRN)
               PYLEFT = PYLAST - PSTRCK * CYRN (IRN)
               PZLEFT = PZLAST - PSTRCK * CZRN (IRN)
               PPLEF2 = PXLEFT**2 + PYLEFT**2 + PZLEFT**2
               UMO2   = ELEFT**2 - PPLEF2
               DELTU2 = AMNRE2 - UMO2
*  |  |  |  +----------------------------------------------------------*
*  |  |  |  |
               IF ( DELTU2 .GT. 0.D+00 ) THEN
*  |  |  |  |  +-------------------------------------------------------*
*  |  |  |  |  |
                  IF ( IUMO .EQ. 0 .AND. PSTRCK .GE. PPLAST ) THEN
                     CXRN (IRN) = CXXINC
                     CYRN (IRN) = CYYINC
                     CZRN (IRN) = CZZINC
                     IUMO   = IUMO + 1
                     GO TO 2020
                  END IF
*  |  |  |  |  |
*  |  |  |  |  +-------------------------------------------------------*
                  PLDOTU = PXLEFT * CXRN (IRN) + PYLEFT * CYRN (IRN) +
     &                     PZLEFT * CZRN (IRN)
                  DELTAE = 0.51D+00 * DELTU2 / ( ELEFT - ( T2 + TSTRCK )
     &                   * PLDOTU / PSTRCK )
*  |  |  |  |  +-------------------------------------------------------*
*  |  |  |  |  |
                  IF ( DELTAE .GE. TSTRCK .OR. IUMO .GT. 3 ) THEN
                     IRN   = IRN - 1
                     KNREJE = KNREJE + 1
                     ANOW  = ANOW + 1.D+00
*  |  |  |  |  |  +----------------------------------------------------*
*  |  |  |  |  |  |
                     IF ( KNREJE .GT. 3 ) THEN
                        LNCASC = .FALSE.
*  |  |  |  |  |  |  +-------------------------------------------------*
*  |  |  |  |  |  |  |
                        IF ( LPCASC ) THEN
                           TPK = TPK + TNK + TN
*  |  |  |  |  |  |  |
*  |  |  |  |  |  |  +-------------------------------------------------*
*  |  |  |  |  |  |  |
                        ELSE
                           TV = TV + TNK + TN
                        END IF
*  |  |  |  |  |  |  |
*  |  |  |  |  |  |  +-------------------------------------------------*
                        TNK = 0.D+00
*  |  |  |  |  |  |
*  |  |  |  |  |  +----------------------------------------------------*
*  |  |  |  |  |  |
                     ELSE
                        TNK = TNK + TN
                     END IF
*  |  |  |  |  |  |
*  |  |  |  |  |  +----------------------------------------------------*
                     GO TO 1900
*  |  |  |  |  |
*  +-<|-<|-<|-<|--<--<--<--<--<
*  |  |  |  |  |
*  |  |  |  |  +-------------------------------------------------------*
*  |  |  |  |  |
                  ELSE
                     TSTRCK = TSTRCK - DELTAE
                     PSTRCK = SQRT ( TSTRCK * ( TSTRCK + 2.D+00 * T2 ) )
                     TKRECL = TKRECL + DELTAE
                     ELEFT  = ELEFT  + DELTAE
                     GO TO 2020
*  |  |  |  |
*  |  |  +-<|-<|--<--<--<--<
                  END IF
*  |  |  |  |  |
*  |  |  |  |  +-------------------------------------------------------*
               END IF
*  |  |  |  |
*  |  |  |  +----------------------------------------------------------*
*  |  |  |
*  |  |  +-------------------------------------------------------------*
*  |  |
*  |  +----------------------------------------------------------------*
*  |  |
         ELSE
            CALL TTRANS ( CX, CY, CZ, COD, SID, SFE, CFE, CXRN (IRN),
     &                    CYRN (IRN), CZRN (IRN) )
         END IF
*  |  |
*  |  +----------------------------------------------------------------*
         IGREYN = IGREYN + 1
         EINCN  = EINCN + TKIN
         ITRN(IRN) = 8
*
*  T2 is the mass energy of the neutron (see before)
*
         TVGRE0 = TVGRE0
         TVGREY = TVGREY + TKIN + TKRECL - TSTRCK - EBNDNG (2)
         EEX = TKIN + TKRECL - TSTRCK - EBNDNG (2)
         TV  = TV + EEX
         ELR (IRN)  = TSTRCK + T2
         PLR(IRN) = PSTRCK
*  |   updating the conservation counters in common balanc
         EINTR  = EINTR + ELR(IRN)
         PXINTR = PXINTR + PLR(IRN)*CXRN(IRN)
         PYINTR = PYINTR + PLR(IRN)*CYRN(IRN)
         PZINTR = PZINTR + PLR(IRN)*CZRN(IRN)
         ICINTR = ICINTR + ICH(ITRN(IRN))
         IBINTR = IBINTR + IBAR(ITRN(IRN))
         IRNTOT = IRNTOT + 1
      GO TO 1900
*  |
*  +--<--<--<--<--< Try to select another cascade nucleon!!!!
 
*  +-------------------------------------------------------------------*
*  |    Cascade protons selection!!!!!
   30 CONTINUE
         IF ( ZNOW .LT. 0.5D+00 .OR. TPK .LE. TPKTE ) THEN
            LPCASC = .FALSE.
            TV  = TV + TPK
            TPK = 0.D+00
            GO TO 1900
         END IF
C
C--------------------------------------------------
C
C   SIMULATION OF CASCADE PROTONS
C
C--------------------------------------------------
C*** PROTON NUMBER  AND PROTON CASCADE ENERGY TEST
C
C
C--------------------------------------------------
C*** SAMPLING OF KINETIC ENERGY TN AND TETA ANGLE DN OF THE EMITTED P
C--------------------------------------------------
C--------------------------------------------------
         T2 = AM(1)
         CALL RBKEKV ( 1, EXSOP, TOEFF, BBTAR, TKIN, TSTRCK,
     &                 PSTRCK, ANOW, TKRECL, COD, SID )
         TN = TKIN
         TPK = TPK - TN
C
C--------------------------------------------------
C*** TEST: IF THE REMAINING ENERGY AND CORRESP.NUCLEON NUMBER ALLOWS FOR
C*** THE CORRESP. NUCLEON EMISSION, ELSE GIVE THE ENERGY OR TAKE
C*** THE MISSING FROM THE EXCITATION
C*** IF FOR LAST CHOSEN NUCLEON TN>REMAINING REST, TAKE IN 90% THE ENERG
C*** DIFFERENCE FROM EXCIT.ENERGY, IN 10% OR IF TV WOULD BE LESS THEN 0.
C***  ADD TN TO TV
C--------------------------------------------------
         IF ( ZNOW .LT. 1.5D+00 .OR. TPK .LE. ETHR  ) THEN
            LPCASC = .FALSE.
            CALL GRNDM(RNDM,1)
            VT  = RNDM(1)
            THK = TV + TPK
            IF ( VT .GT. 0.9D+00 .OR. THK .LT. 0.D+00 ) THEN
               TV=TV+TPK+TN
               TPK = 0.D+00
               GO TO 1900
            END IF
            TPK = 0.D0
            TV  = THK
         END IF
C
C--------------------------------------------------
C*** CORRESP. NUCLEON NUMBER COUNTING
C--------------------------------------------------
         IF (TN .LT. ETHR) THEN
            TV = TV + TN
            GO TO 1900
*  |  |
*  +-<|--<--<--<--<--< try to select another nucleon
         END IF
         ANOW = ANOW - 1.D0
         ZNOW = ZNOW - 1.D0
         FRSURV = ANOW / BBTAR
C
C--------------------------------------------------
C*** TAKE THE EMITTED CORRESP. NUCLEON INTO THE FINAL PARTICLE TABLE
C*** C./FINUC/ AFTER CHOICE OF PHI AND TURNING INTO THE LAB SYSTEM
C--------------------------------------------------
         IRN = IRN + 1
         CALL COSI ( SFE, CFE )
*  |  +----------------------------------------------------------------*
*  |  |
         IF ( ANOW .LT. 30.D+00 .OR. FRSURV .LT. 0.66D+00 ) THEN
            PXLAST = PXTTOT - PXINTR - PXNUCR
            PYLAST = PYTTOT - PYINTR - PYNUCR
            PZLAST = PZTTOT - PZINTR - PZNUCR
            PPLAS2 = PXLAST**2 + PYLAST**2 + PZLAST**2
*  |  |  +-------------------------------------------------------------*
*  |  |  |
            IF ( PPLAS2 .GT. 1.D-6 * PTTOT**2 ) THEN
               PPLAST = SQRT ( PPLAS2 )
               CXXINC = PXLAST / PPLAST
               CYYINC = PYLAST / PPLAST
               CZZINC = PZLAST / PPLAST
               CDMIN = MIN ( ONEONE, PSTRCK / PPLAST ) * ( 1.D+00 -
     &                       0.8D+00 * FRSURV )
               COD = 1.D+00 - 0.5D+00 * ( 1.D+00 - COD ) * ( 1.D+00 -
     &               CDMIN )
               COD = MIN ( COD, ONEONE )
               SID = SQRT ( 1.D+00 - COD**2 )
*  |  |  |
*  |  |  +-------------------------------------------------------------*
*  |  |  |
            ELSE
               CXXINC = CX
               CYYINC = CY
               CZZINC = CZ
            END IF
*  |  |  |
*  |  |  +-------------------------------------------------------------*
            CALL TTRANS ( CXXINC, CYYINC, CZZINC, COD, SID, SFE, CFE,
     &                    CXRN (IRN), CYRN (IRN), CZRN (IRN) )
*  |  |  +-------------------------------------------------------------*
*  |  |  |
            IF ( NINT (ANOW) .GT. 0 ) THEN
               AMNRES = AMUAMU * ANOW - ZNOW * AMELEC + 1.D-03 * FKENER
     &                ( ANOW, ZNOW ) + ELBNDE ( NINT ( ZNOW ) )
*  |  |  |
*  |  |  +-------------------------------------------------------------*
*  |  |  |
            ELSE
               AMNRES = 0.D+00
            END IF
*  |  |  |
*  |  |  +-------------------------------------------------------------*
            AMNRE2 = AMNRES  * AMNRES
            ELEFT  = ETTOT  - ENUCR  - EINTR  - TSTRCK - T2
            IUMO   = 0
*  |  |  +-------------------------------------------------------------*
*  |  |  |
3030        CONTINUE
               PXLEFT = PXLAST - PSTRCK * CXRN (IRN)
               PYLEFT = PYLAST - PSTRCK * CYRN (IRN)
               PZLEFT = PZLAST - PSTRCK * CZRN (IRN)
               PPLEF2 = PXLEFT**2 + PYLEFT**2 + PZLEFT**2
               UMO2   = ELEFT**2 - PPLEF2
               DELTU2 = AMNRE2 - UMO2
*  |  |  |  +----------------------------------------------------------*
*  |  |  |  |
               IF ( DELTU2 .GT. 0.D+00 ) THEN
*  |  |  |  |  +-------------------------------------------------------*
*  |  |  |  |  |
                  IF ( IUMO .EQ. 0 .AND. PSTRCK .GE. PPLAST ) THEN
                     CXRN (IRN) = CXXINC
                     CYRN (IRN) = CYYINC
                     CZRN (IRN) = CZZINC
                     IUMO = IUMO + 1
                     GO TO 3030
                  END IF
*  |  |  |  |  |
*  |  |  |  |  +-------------------------------------------------------*
                  IUMO   = IUMO + 1
                  PLDOTU = PXLEFT * CXRN (IRN) + PYLEFT * CYRN (IRN) +
     &                     PZLEFT * CZRN (IRN)
                  DELTAE = 0.51D+00 * DELTU2 / ( ELEFT - ( T2
     &                   + TSTRCK ) * PLDOTU / PSTRCK )
*  |  |  |  |  +-------------------------------------------------------*
*  |  |  |  |  |
                  IF ( DELTAE .GE. TSTRCK .OR. IUMO .GT. 3 ) THEN
                     IRN   = IRN - 1
                     KPREJE = KPREJE + 1
                     ANOW  = ANOW + 1.D+00
                     ZNOW  = ZNOW + 1.D+00
*  |  |  |  |  |  +----------------------------------------------------*
*  |  |  |  |  |  |
                     IF ( KPREJE .GT. 3 ) THEN
                        LPCASC = .FALSE.
*  |  |  |  |  |  |  +-------------------------------------------------*
*  |  |  |  |  |  |  |
                        IF ( LNCASC ) THEN
                           TNK = TPK + TNK + TN
*  |  |  |  |  |  |  |
*  |  |  |  |  |  |  +-------------------------------------------------*
*  |  |  |  |  |  |  |
                        ELSE
                           TV = TV + TPK + TN
                        END IF
*  |  |  |  |  |  |  |
*  |  |  |  |  |  |  +-------------------------------------------------*
                        TPK = 0.D+00
*  |  |  |  |  |  |
*  |  |  |  |  |  +----------------------------------------------------*
*  |  |  |  |  |  |
                     ELSE
                        TPK = TPK + TN
                     END IF
*  |  |  |  |  |  |
*  |  |  |  |  |  +----------------------------------------------------*
                     GO TO 1900
*  |  |  |  |  |
*  +-<|-<|-<|-<|--<--<--<--<--<
*  |  |  |  |  |
*  |  |  |  |  +-------------------------------------------------------*
*  |  |  |  |  |
                  ELSE
                     TSTRCK = TSTRCK - DELTAE
                     PSTRCK = SQRT ( TSTRCK * ( TSTRCK + 2.D+00 * T2 ) )
                     TKRECL = TKRECL + DELTAE
                     ELEFT  = ELEFT  + DELTAE
                     GO TO 3030
*  |  |  |  |  |
*  |  |  +-<|-<|--<--<--<--<
                  END IF
*  |  |  |  |  |
*  |  |  |  |  +-------------------------------------------------------*
               END IF
*  |  |  |  |
*  |  |  |  +----------------------------------------------------------*
*  |  |  |
*  |  |  +-------------------------------------------------------------*
*  |  |
*  |  +----------------------------------------------------------------*
*  |  |
         ELSE
            CALL TTRANS ( CX, CY, CZ, COD, SID, SFE, CFE, CXRN (IRN),
     &                    CYRN (IRN), CZRN (IRN) )
         END IF
*  |  |
*  |  +----------------------------------------------------------------*
         IGREYP = IGREYP + 1
         EINCP  = EINCP + TKIN
C
         ITRN (IRN) = 1
*
*  T2 is the mass energy of the proton (see before)
*
*        TVGRE0 = TVGRE0 + EXSOP (1)
         TVGRE0 = TVGRE0
         TVGREY = TVGREY + TKIN + TKRECL - TSTRCK - EBNDNG (1)
         EEX = TKIN + TKRECL - TSTRCK - EBNDNG (1)
         TV  = TV + EEX
         ELR (IRN) = TSTRCK + T2
         PLR (IRN) = PSTRCK
*  |   updating the conservation counters in common balanc
         EINTR  = EINTR + ELR(IRN)
         PXINTR = PXINTR + PLR(IRN)*CXRN(IRN)
         PYINTR = PYINTR + PLR(IRN)*CYRN(IRN)
         PZINTR = PZINTR + PLR(IRN)*CZRN(IRN)
         ICINTR = ICINTR + ICH(ITRN(IRN))
         IBINTR = IBINTR + IBAR(ITRN(IRN))
         IRNTOT = IRNTOT + 1
      GO TO 1900
*  |
*  +--<--<--<--<--< Try to select another cascade nucleon!!!!
C
C--------------------------------------------------
C*** END OF CASCADE NUCLEON EMISSION
C*** HANDLING OF THE CASES: NO OR ONE PARTICLE IN H-A-FINAL STATE
C*** FINAL CALCULATION OF TV
C--------------------------------------------------
2000  CONTINUE
*  +-------------------------------------------------------------------*
*  |   Now working with atomic masses !
      IF ( NINT (ANOW) .GT. 0 ) THEN
         AMMRES = AMUAMU * ANOW + 1.D-03 * FKENER ( ANOW, ZNOW )
         AMNRES = AMMRES - ZNOW * AMELEC + ELBNDE ( NINT (ZNOW) )
         DELEFT = AMMTAR - AMNTAR - ( ICHTAR - ZNOW ) * AMELEC
*  |
*  +-------------------------------------------------------------------*
*  |
      ELSE
         AMMRES = 0.D+00
         AMNRES = 0.D+00
         DELEFT = 0.D+00
      END IF
*  |
*  +-------------------------------------------------------------------*
*  +-------------------------------------------------------------------*
*  |                        Check the number of emitted particles
      IF ( IRN - NP0 - 1 ) 4001,4002,4003
4001  CONTINUE
*  |  No particle emitted !!!
*  |  Now working with atomic masses
         TV = ETTOT - AMNRES
      GO TO 4013
*  |
*  +-------------------------------------------------------------------*
*  |  One particle emitted
4002  CONTINUE
*  |  +----------------------------------------------------------------*
*  |  |
         IF ( TV .LE. 0.001D+00 ) THEN
            TV  = TVEUZ0
         END IF
*  |  |
*  |  +----------------------------------------------------------------*
*  |   Now working with atomic masses
         UMO2   = ( ETTOT + DELEFT )**2 - PTTOT**2
         IUMO   = 0
*  |  +----------------------------------------------------------------*
*  |  |
4040     CONTINUE
            AMSTAR = AMMRES + TV
            AMEMIT = AM ( ITRN (NP0+1) )
            UMIN2  = ( AMEMIT + AMSTAR )**2
*  |  |  +-------------------------------------------------------------*
*  |  |  |
            IF ( UMO2 .LE. UMIN2 ) THEN
               DELTAE = 0.51D+00 * ( UMIN2 - UMO2 ) / ( AMEMIT +
     &                  AMSTAR )
               IUMO = IUMO + 1
*  |  |  |  +----------------------------------------------------------*
*  |  |  |  |
               IF ( TV .GE. DELTAE .AND. IUMO .LE. 3 ) THEN
                  TV = TV - DELTAE
                  GO TO 4040
*  |  |  |  |
*  |  +-<|-<|--<--<--<--<--<
*  |  |  |  |
*  |  |  |  +----------------------------------------------------------*
*  |  |  |  |
               ELSE
                  WRITE ( LUNERR, * ) ' Nucriv: single particle emiss',
     &                  'ion not possible for missing invariant mass!!',
     &                     NINT (ANOW), NINT (ZNOW), UMIN2, UMO2,
     &                     ITRN (NP0+1), IUMO
                  ANOW = BBTAR
                  ZNOW = ZZTAR
                  KTARP = 0
                  KTARN = 0
                  IGREYP = 0
                  IGREYN = 0
                  ENUCR  = 0.D+00
                  PXNUCR = 0.D+00
                  PYNUCR = 0.D+00
                  PZNUCR = 0.D+00
                  EINTR   = 0.D+00
                  PXINTR  = 0.D+00
                  PYINTR  = 0.D+00
                  PZINTR  = 0.D+00
                  TVGREY  = 0.D+00
                  UMO     = SQRT ( UMO2 )
                  TVGRE0  = UMO - AMMTAR
                  TV  = TVGRE0
                  IRN = NP0
                  RETURN
               END IF
*  |  |  |  |
*  |  |  |  +----------------------------------------------------------*
            END IF
*  |  |  |
*  |  |  +-------------------------------------------------------------*
*  |  |
*  |  +----------------------------------------------------------------*
         UMO    = SQRT ( UMO2 )
         GAMCM  = ( ETTOT + DELEFT ) / UMO
         ETACM  = PTTOT / UMO
         ETAX   = PXTTOT / UMO
         ETAY   = PYTTOT / UMO
         ETAZ   = PZTTOT / UMO
         PXNUCR = PLR (NP0+1) * CXRN (NP0+1)
         PYNUCR = PLR (NP0+1) * CYRN (NP0+1)
         PZNUCR = PLR (NP0+1) * CZRN (NP0+1)
         ETAPCM = - ETAX * PXNUCR - ETAY * PYNUCR - ETAZ * PZNUCR
         PHELP  = ETAPCM / ( GAMCM + 1.D+00 ) + ELR (NP0+1)
*  |  Old direction cosines of the emitted particle in CMS
         PXINTR = PXNUCR - ETAX * PHELP
         PYINTR = PYNUCR - ETAY * PHELP
         PZINTR = PZNUCR - ETAZ * PHELP
         PCMS   = SQRT ( PXINTR**2 + PYINTR**2 + PZINTR**2 )
         CXXINC = PXINTR / PCMS
         CYYINC = PYINTR / PCMS
         CZZINC = PZINTR / PCMS
         ECMSST = 0.5D+00 * ( UMO2 + AMSTAR**2 - AMEMIT**2 ) / UMO
         ECMSEM = UMO - ECMSST
         PCMS   = SQRT ( ECMSEM**2 - AMEMIT**2 )
*  |   Now save the old direction
         PXINTR = PCMS * CXXINC
         PYINTR = PCMS * CYYINC
         PZINTR = PCMS * CZZINC
         ETAPCM = ETAX * PXINTR + ETAY * PYINTR + ETAZ * PZINTR
         ELR (NP0+1) = GAMCM * ECMSEM + ETAPCM
         PHELP   = ETAPCM / ( GAMCM + 1.D+00 ) + ECMSEM
         PXNUCR  = PXINTR + ETAX * PHELP
         PYNUCR  = PYINTR + ETAY * PHELP
         PZNUCR  = PZINTR + ETAZ * PHELP
         PLR (NP0+1) = SQRT ( PXNUCR**2 + PYNUCR**2 + PZNUCR**2 )
         ENUCR   = ELR  (NP0+1)
         CXRN (NP0+1) = PXNUCR / PLR (NP0+1)
         CYRN (NP0+1) = PYNUCR / PLR (NP0+1)
         CZRN (NP0+1) = PZNUCR / PLR (NP0+1)
         KTARP   = KTARP + IGREYP
         KTARN   = KTARN + IGREYN
         IGREYP  = 0
         IGREYN  = 0
         EINTR   = 0.D+00
         PXINTR  = 0.D+00
         PYINTR  = 0.D+00
         PZINTR  = 0.D+00
         TVGRE0  = 0.D+00
         TVGREY  = 0.D+00
         TVEUZ   = TV
      GO TO 4013
*  |
*  +-------------------------------------------------------------------*
 
*  +-------------------------------------------------------------------*
*  |  More than one particle emitted!!
4003  CONTINUE
*  |  Now working with atomic masses !
         TV = ETTOT + DELEFT - AMMRES - ENUCR - EINTR
*  |  +----------------------------------------------------------------*
*  |  |
         IF ( TV .LE. 0.D+00 ) THEN
            WRITE ( LUNERR, * )' *** Nucrin failure: Tv < 0!!',TV
         END IF
*  |  |
*  |  +----------------------------------------------------------------*
 4013 CONTINUE
*  |
*  +-------------------------------------------------------------------*
*  +-------------------------------------------------------------------*
*  |  Check if a neutron/proton has been left: if so add to the
*  |  emitted particle list, conserving energyu (but not momentum)
      IF ( NINT (ANOW) .EQ. 1 ) THEN
         IF ( TV .LT. 0.D+00 ) THEN
            LRESMP = .TRUE.
            RETURN
         END IF
         PXRES = PXTTOT - PXINTR - PXNUCR
         PYRES = PYTTOT - PYINTR - PYNUCR
         PZRES = PZTTOT - PZINTR - PZNUCR
         PTRES = SQRT ( PXRES**2 + PYRES**2 + PZRES**2 )
         IRN = IRN + 1
         ANOW  = ANOW  - 1.D+00
         IF ( NINT (ZNOW) .GT. 0.D+00 ) THEN
            ITRN (IRN) = 1
            ZNOW  = ZNOW  - 1.D+00
            KTARP = KTARP + 1
            IGREYP = IGREYP + 1
            TSTRCK = ETTOT - AM (ITRN(IRN)) - ENUCR - EINTR
            EINCP  = EINCP + TSTRCK
         ELSE
            ITRN (IRN) = 1
            KTARN = KTARN + 1
            IGREYN = IGREYN + 1
            TSTRCK = ETTOT - AM (ITRN(IRN)) - ENUCR - EINTR
            EINCN  = EINCN + TSTRCK
         END IF
         TV  = 0.D+00
         ELR (IRN) = TSTRCK + AM (ITRN(IRN))
         PLR (IRN) = SQRT ( TSTRCK * ( TSTRCK + 2.D+00 * AM(ITRN(IRN)) )
     &             )
         PTRES = PLR (IRN) / MAX ( PTRES, ANGLGB )
         CXRN (IRN) = PXRES * PTRES
         CYRN (IRN) = PYRES * PTRES
         CZRN (IRN) = PZRES * PTRES
*  |   updating the conservation counters in common balanc
         EINTR  = EINTR + ELR(IRN)
         PXINTR = PXINTR + PLR(IRN)*CXRN(IRN)
         PYINTR = PYINTR + PLR(IRN)*CYRN(IRN)
         PZINTR = PZINTR + PLR(IRN)*CZRN(IRN)
         ICINTR = ICINTR + ICH(ITRN(IRN))
         IBINTR = IBINTR + IBAR(ITRN(IRN))
         IRNTOT = IRNTOT + 1
      END IF
*  |
*  +-------------------------------------------------------------------*
      RETURN
      END
