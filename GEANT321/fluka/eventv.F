*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:19:55  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.42  by  S.Giani
*-- Author :
*$ CREATE EVENTV.FOR
*COPY EVENTV
*
*=== eventv ===========================================================*
*
      SUBROUTINE EVENTV ( IJIJ, POO, EKE, TXI, TYI, TZI, WE, IMAT )
 
#include "geant321/dblprc.inc"
#include "geant321/dimpar.inc"
#include "geant321/iounit.inc"
*
*----------------------------------------------------------------------*
*                                                                      *
*  Created  on  7 september 1989      by         Alfredo Ferrari       *
*                                                 INFN - Milan         *
*                                                                      *
*  Last change  on  15 april 1993     by         Alfredo Ferrari       *
*                                                                      *
*  This routine is a strongly modified and improved version of the     *
*  original Eventq routine of the Fluka package: it requires modified  *
*  versions of many other routines                                     *
*                                                                      *
*  Eventv90: new version by A. Ferrari, INFN-Milan and CERN-SPS, 7-9-89*
*       main modifications:                                            *
*                  - new treatment for hydrogen (following Moehring's  *
*                                                suggestion)           *
*                  - new commons and coding for energy, momentum,      *
*                    electric and baryonic charge conservation         *
*                  - new correlations both in the "low" (p < 5 GeV/c)  *
*                    and in the "high" energy models                   *
*                  - complete revision of kinematics for both models   *
*                  - tentative treatment for Lambdas, Sigmas with the  *
*                    "high" energy model down to 2.5 GeV/c to overcome *
*                    the impossibility of Hadrin to treat these parti- *
*                    cles                                              *
*                  - tentative introduction of Xsi0, Xsi-, Omega and   *
*                    their antiparticles and of Sigmabars's            *
*                  - complete new treatment of cascade nucleons        *
*                  - accurate treatment of binding energy effects,     *
*                    using atomic mass tables and the proper isotopic  *
*                    composition of elements                           *
*                  - possibility to use an evaporation module (derived *
*                    from the Hetc-Kfa one) after the interaction      *
*                  - possibility to use a nuclear gamma deexcitation   *
*                    module (written by P.Sala, INFN-Milan) to produce *
*                    deexcitation gammas after the evaporation step    *
*                                                                      *
*----------------------------------------------------------------------*
*
      PARAMETER ( COS45 = 0.7071067811865475D+00 )
      PARAMETER ( SIN30 = 0.5D+00 )
      PARAMETER ( COS30 = 0.8660254037844387D+00 )
      PARAMETER ( CSNPRN = 50.D+00 * CSNNRM )
#include "geant321/balanc.inc"
#include "geant321/corinc.inc"
#include "geant321/depnuc.inc"
#include "geant321/fheavy.inc"
#include "geant321/finlsp3.inc"
#include "geant321/finuc.inc"
#include "geant321/hadflg.inc"
#include "geant321/hadpar.inc"
#include "geant321/isotop.inc"
#include "geant321/mapa.inc"
#include "geant321/nucdat.inc"
#include "geant321/nucpar.inc"
#include "geant321/part.inc"
#include "geant321/parevt.inc"
#include "geant321/resnuc.inc"
      COMMON / FKCASF / PKFRMI, COSTH, PKIN
      COMMON /FKEVNT/ LNUCRI, LHADRI
      LOGICAL LNUCRI, LHADRI, LINCCK, LACCEP, LCHTYP
      REAL RNDM(1)
*  Note Pthrsh is also in Ferevv!!
      DIMENSION SOPPP(2),EXSOP(2),PTHRSH(39),IJNUCR(39)
* Use this "save" if common dbgdbg is commented
      SAVE PTHRSH, IPRI, IEVT, IJNUCR
      DATA PTHRSH / 16*5.D+00,2*2.5D+00,5.D+00,3*2.5D+00,8*5.D+00,
     &              9*2.5D+00 /
      DATA IJNUCR / 16*1,2*0,1,3*0,8*1,9*0 /
      DATA IPRI / 0 /
      DATA IEVT / 0 /
*
*  +-------------------------------------------------------------------*
*  |  Heavy ions are treated in eventa: now switched off!!
      IF (IJIJ .EQ. 30) THEN
         STOP 'HEAVY ION'
      END IF
*  |
*  +-------------------------------------------------------------------*
*   First of all: get the "part" index of the incoming "paprop" ijij
      IJ = IPTOKP (IJIJ)
      IEVT = IEVT + 1
      AMCHCK = 0.5D+00 * ( POO * POO - EKE * EKE ) / EKE
      AHELP  =  ABS ( AMCHCK - AM ( IJ ) )
*  +-------------------------------------------------------------------*
*  |
      IF ( AHELP .GT. ANGLGB * AM ( IJ ) ) THEN
         POO = SQRT ( EKE * ( EKE + 2.D+00 * AM (IJ) ) )
      END IF
*  |
*  +-------------------------------------------------------------------*
*
*  Variable initialization for conservation checks
*
      PTTOT  = POO
      PXTTOT = PTTOT*TXI
      PYTTOT = PTTOT*TYI
      PZTTOT = PTTOT*TZI
      ICHTAR = NINT(ZTAR(IMAT))
      ICHTOT = ICHTAR + ICH(IJ)
*  +-------------------------------------------------------------------*
*  |    Choice of the mass number of the target nucleus: use the input
*  |    one if any
      IF ( MSSNUM (IMAT) .GT. 0 ) THEN
         IBTAR = MSSNUM (IMAT)
*  |
*  +-------------------------------------------------------------------*
*  |    Choice of the mass number of the target nucleus according
*  |    to the natural isotopic composition of element ichtar
      ELSE
         CALL GRNDM(RNDM,1)
         RNDISO = RNDM (1)
*  |  +----------------------------------------------------------------*
*  |  |    Loop on the stable isotopes
         DO 25 IS = ISONDX (1,ICHTAR), ISONDX (2,ICHTAR) - 1
            RNDISO = RNDISO - ABUISO (IS)
            IF ( RNDISO .LT. 0.D+00 ) GO TO 30
   25    CONTINUE
*  |  |
*  |  +----------------------------------------------------------------*
         IS = ISONDX (2,ICHTAR)
   30    CONTINUE
         IBTAR = ISOMNM (IS)
      END IF
*  |
*  +-------------------------------------------------------------------*
      IBTOT  = IBTAR + IBAR(IJ)
      BBTAR  = IBTAR
      ZZTAR  = ICHTAR
      ZTO103 = ZZTAR**0.3333333333333333D+00
*  +-------------------------------------------------------------------*
*  |
      IF ( IBTAR .EQ. 1 ) THEN
         ITTA   = 8 - 7 * ICHTAR
         ETTOT  = AM (ITTA) + EKE + AM(IJ)
         AMNTAR = AM (ITTA)
*  |
*  +-------------------------------------------------------------------*
*  | The following should be done with the proper mass of the nuclide
      ELSE
*        AMNTAR = BBTAR * AMUC12
         AMMTAR = BBTAR * AMUAMU + 0.001D+00 * FKENER ( BBTAR, ZZTAR )
         AMNTAR = AMMTAR - ZZTAR * AMELEC + ELBNDE ( ICHTAR )
         ETTOT  = AMNTAR + EKE + AM (IJ)
      END IF
*  |
*  +-------------------------------------------------------------------*
*  +-------------------------------------------------------------------*
*  |  Kaon-short and kaon-long are transformed into a kaon0 or a
*  |  kaon0-bar
      IF (IJ .EQ. 12 .OR. IJ. EQ. 19) THEN
         IJ = 24
         CALL GRNDM(RNDM,1)
         IF (RNDM(1) .GT. 0.5D0) IJ = 25
*  |
*  +-------------------------------------------------------------------*
*  |  Pi0 quark configurations are selected according to 50% uubar and
*  |  50% ddbar
      ELSE IF ( IJ .EQ. 23 .OR. IJ .EQ. 26 ) THEN
*  |  +----------------------------------------------------------------*
*  |  |
         CALL GRNDM(RNDM,1)
         IF ( RNDM (1) .LT. 0.5D+00 ) THEN
            IJ = 23
*  |  |
*  |  +----------------------------------------------------------------*
*  |  |
         ELSE
            IJ = 26
         END IF
*  |  |
*  |  +----------------------------------------------------------------*
      END IF
*  |
*  +-------------------------------------------------------------------*
      ETEPS  = MAX ( 1.D-10 * EKE, 1.D-10 )
  50  CONTINUE
*-->-->-->-->--> Come here for resampling
* Reset all the accumulators
      NP     = NP0
      NPHEAV = 0
      LRNFSS = .FALSE.
      LEVDIF = .FALSE.
      LHADRI = .FALSE.
      ICHTOT = ICHTAR + ICH(IJ)
      TV = 0.D+00
      TVRECL = 0.D+00
      TVHEAV = 0.D+00
      TVCMS  = 0.D+00
      EOTEST = ETTOT
      LRESMP = .FALSE.
      LLASTN = .FALSE.
      LLAST1 = .FALSE.
      EUZ = 0.D+00
      PUX = 0.D+00
      PUY = 0.D+00
      PUZ = 0.D+00
      TVEUZ = 0.D+00
      ICU = 0
      IBU = 0
      EINTR  = 0.D+00
      PXINTR = 0.D+00
      PYINTR = 0.D+00
      PZINTR = 0.D+00
      ICINTR = 0
      IBINTR = 0
      ENUCR  = 0.D+00
      PXNUCR = 0.D+00
      PYNUCR = 0.D+00
      PZNUCR = 0.D+00
      ICNUCR = 0
      IBNUCR = 0
      EFRM  = 0.D+00
      PXFRM = 0.D+00
      PYFRM = 0.D+00
      PZFRM = 0.D+00
      PSEA  = 0.D+00
      TVGRE0  = 0.D+00
      TVGREY  = 0.D+00
      IGREYP  = 0
      IGREYN  = 0
      KTARP = 0
      KTARN = 0
      IEVAPL = 0
      IEVAPH = 0
      IEVPRO = 0
      IEVNEU = 0
      IEVDEU = 0
      IEVTRI = 0
      IEV3HE = 0
      IEV4HE = 0
      IDEEXG = 0
      ANOW   = BBTAR
      ZNOW   = ZZTAR
      DSOPP  = 0.D+00
      NNHAD  = 0
      RN1GSC = -1.D+00
      RN2GSC = -1.D+00
      IF ( POO .GE. PTHRSH (IJIJ) ) GO TO 200
* Below 5 GeV/c use Nucrin:
 100  CONTINUE
      ANCOLL = 1.D+00
*  +-------------------------------------------------------------------*
*  |  Check for Hydrogen
      IF ( IBTAR .NE. 1  ) THEN
*  |  +----------------------------------------------------------------*
*  |  |  Steering for Peanut: very temporary
*  |  |  Protons and Neutrons below 260 MeV:
         IF ( EKE .LT. 0.260D+00 .AND. ( IJ .EQ. 1 .OR. IJ .EQ. 8 )
     &        .AND. IBTAR .GT. 4 ) THEN
            CALL PEANUT ( IJ, EKE, POO, TXI, TYI, TZI, WE )
            IF ( LRESMP ) GO TO 50
            RETURN
*  |  |
*  |  +----------------------------------------------------------------*
*  |  |  Stopping pi-s:
         ELSE IF ( IJ .EQ. 14 .AND. EKE .LT. 2.D+00 * GAMMIN .AND. IBTAR
     &             .GT. 4 ) THEN
            CALL PEANUT ( IJ, EKE, POO, TXI, TYI, TZI, WE )
            IF ( LRESMP ) GO TO 50
            RETURN
         END IF
*  |  |
*  |  +----------------------------------------------------------------*
*  |  Now at first we sample the number of interactions
         CALL CORRIN ( ZZTAR, BBTAR, IJ, POO, EKE )
      END IF
      LNUCRI = .TRUE.
      NP = NP0
      ELAB  = EKE + AM(IJ)
*  Redefinition of particle types: the recovery of proper charge and
*  strangeness conservation after the interaction is not yet implemented
      KPROJ = IJ
      LCHTYP = .FALSE.
*  +-------------------------------------------------------------------*
*  |
      IF ( KPROJ .LE. 30 ) THEN
         GO TO ( 2017, 2018, 2999, 2020, 2021, 2022, 2999, 2999, 2999,
     &           2026 ), KPROJ - 16
         GO TO 2999
*  |
*  +-------------------------------------------------------------------*
*  |
      ELSE
         GO TO ( 2031, 2032, 2033, 2034, 2035, 2036, 2037, 2038, 2039 ),
     &      KPTOIP (KPROJ) - 30
         GO TO 2999
      END IF
*  |
*  +-------------------------------------------------------------------*
*  +-------------------------------------------------------------------*
*  |  Lambdas are considered as neutrons: after the interaction
*  |  a possible way to recover strangeness and charge conservation
*  |  is to flip a pi- --> K-, or an eventual neutron into
*  |  a Sigma0/Lambda,
*  |  (we must flip a d quark into an s quark), or a pi0 --> K0bar,
*  |  or a p --> Sigma+
 2017 CONTINUE
         LCHTYP = .TRUE.
         KPROJ = 8
*  |   Adjusting mass and momentum for lambdas (now n)
         POO = SQRT ( ( ELAB + AM(KPROJ) ) * ( ELAB - AM(KPROJ) ) )
         PTTOT  = POO
         PXTTOT = PTTOT*TXI
         PYTTOT = PTTOT*TYI
         PZTTOT = PTTOT*TZI
      GO TO 2999
*  |
*  +-------------------------------------------------------------------*
*  |  Alambdas are considered as antineutrons : after the interaction
*  |  a possible way to recover strangeness and charge conservation
*  |  is to flip a pi+ --> K+, or an eventual neutronbar into a
*  |  Lambdabar, (we must flip a dbar quark into an sbar quark),
*  |  or a pi0 --> K0, or a pbar --> Sigma+bar (ASigma-) ( there
*  |  are also possible ways through eta, rho and omega mesons )
 2018 CONTINUE
         LCHTYP = .TRUE.
         KPROJ = 9
*  |   Adjusting mass and momentum for alambdas (now an)
         POO = SQRT ( ( ELAB + AM(KPROJ) ) * ( ELAB - AM(KPROJ) ) )
         PTTOT  = POO
         PXTTOT = PTTOT*TXI
         PYTTOT = PTTOT*TYI
         PZTTOT = PTTOT*TZI
      GO TO 2999
*  |
*  +-------------------------------------------------------------------*
*  |  Sigma-s are considered as neutrons : after the interaction
*  |  a possible way to recover strangeness and charge conservation
*  |  is to flip a pi+ --> K0bar, or an eventual neutron into a Sigma-,
*  |  (we must flip a u quark into an s quark), or a pi0 --> K-,
*  |  or a p --> Sigma0/Lambda
 2020 CONTINUE
         LCHTYP = .TRUE.
         KPROJ = 8
*  |   The following card must be commented following the strangeness
*  |   and charge conservation recovery implementation
         ICHTOT = ICHTOT + 1
*  |   Adjusting mass and momentum for Sigma-s (now n)
         POO = SQRT ( ( ELAB + AM(KPROJ) ) * ( ELAB - AM(KPROJ) ) )
         PTTOT  = POO
         PXTTOT = PTTOT*TXI
         PYTTOT = PTTOT*TYI
         PZTTOT = PTTOT*TZI
      GO TO 2999
*  |
*  +-------------------------------------------------------------------*
*  |  Sigma+s are considered as protons : after the interaction
*  |  a possible way to recover strangeness and charge conservation
*  |  is to flip a pi- --> K-, or an eventual proton into a Sigma+,
*  |  (we must flip a d quark into an s quark), or a pi0 --> K0bar,
*  |  or a n --> Sigma0/Lambda
 2021 CONTINUE
         LCHTYP = .TRUE.
         KPROJ = 1
*  |   Adjusting mass and momentum for Sigma+s (now p)
         PTTOT  = POO
         PXTTOT = PTTOT*TXI
         PYTTOT = PTTOT*TYI
         PZTTOT = PTTOT*TZI
      GO TO 2999
*  |
*  +-------------------------------------------------------------------*
*  |  Sigma0s are considered as neutrons: after the interaction
*  |  a possible way to recover strangeness and charge conservation
*  |  is to flip a pi- --> K-, or an eventual neutron into
*  |  a Sigma0/Lambda,
*  |  (we must flip a d quark into an s quark), or a pi0 --> K0bar,
*  |  or a p --> Sigma+
 2022 CONTINUE
         LCHTYP = .TRUE.
         KPROJ = 8
*  |   Adjusting mass and momentum for Sigma0s (now n)
         POO = SQRT ( ( ELAB + AM(KPROJ) ) * ( ELAB - AM(KPROJ) ) )
         PTTOT  = POO
         PXTTOT = PTTOT*TXI
         PYTTOT = PTTOT*TYI
         PZTTOT = PTTOT*TZI
      GO TO 2999
*  |
*  +-------------------------------------------------------------------*
*  |  Pi0 must be 23 for Nucrin
 2026 CONTINUE
         KPROJ = 23
      GO TO 2999
*  |
*  +-------------------------------------------------------------------*
*  |  ASigma-s are considered as pbars: after the interaction
*  |  a possible way to recover strangeness and charge conservation
*  |  is to flip a pi+ --> K+, or an eventual pbar into a ASigma-,
*  |  (we must flip a dbar quark into an sbar quark), or a pi0 --> K0,
*  |  or a nbar --> ASigma0/ALambda
 2031 CONTINUE
         LCHTYP = .TRUE.
         KPROJ = 2
*  |   Adjusting mass and momentum for ASigma-s (now pbar)
         POO = SQRT ( ( ELAB + AM(KPROJ) ) * ( ELAB - AM(KPROJ) ) )
         PTTOT  = POO
         PXTTOT = PTTOT*TXI
         PYTTOT = PTTOT*TYI
         PZTTOT = PTTOT*TZI
      GO TO 2999
*  |
*  +-------------------------------------------------------------------*
*  |  ASigma0s are considered as nbar: after the interaction
*  |  a possible way to recover strangeness and charge conservation
*  |  is to flip a pi+ --> K+, or an eventual nbar into
*  |  a ASigma0/ALambda,
*  |  (we must flip a dbar quark into an sbar quark), or a pi0 --> K0,
*  |  or a pbar --> ASigma-
 2032 CONTINUE
         LCHTYP = .TRUE.
         KPROJ = 9
*  |   Adjusting mass and momentum for ASigma0s (now nbar)
         POO = SQRT ( ( ELAB + AM(KPROJ) ) * ( ELAB - AM(KPROJ) ) )
         PTTOT  = POO
         PXTTOT = PTTOT*TXI
         PYTTOT = PTTOT*TYI
         PZTTOT = PTTOT*TZI
      GO TO 2999
*  |
*  +-------------------------------------------------------------------*
*  |  ASigma+s are considered as nbar : after the interaction
*  |  a possible way to recover strangeness and charge conservation
*  |  is to flip a pi- --> K0, or an eventual nbar into a ASigma+,
*  |  (we must flip a ubar quark into an sbar quark), or a pi0 --> K+,
*  |  or a pbar --> ASigma0/ALambda
 2033 CONTINUE
         LCHTYP = .TRUE.
         KPROJ = 9
*  |   The following card must be commented following the strangeness
*  |   and charge conservation recovery implementation
         ICHTOT = ICHTOT - 1
*  |   Adjusting mass and momentum for ASigma+s (now nbar)
         POO = SQRT ( ( ELAB + AM(KPROJ) ) * ( ELAB - AM(KPROJ) ) )
         PTTOT  = POO
         PXTTOT = PTTOT*TXI
         PYTTOT = PTTOT*TYI
         PZTTOT = PTTOT*TZI
      GO TO 2999
*  |
*  +-------------------------------------------------------------------*
*  |  Xsi0s are considered as neutrons: after the interaction
*  |  a possible way to recover strangeness and charge conservation
*  |  is to flip an eventual neutron into a Xsi0,
*  |  or a proton --> Sigma+ and a pi- --> K- (or a pi0 --> K0bar) etc.
*  |  (we must flip two d quarks into s quarks)
 2034 CONTINUE
         LCHTYP = .TRUE.
         KPROJ = 8
*  |   Adjusting mass and momentum for Xsi0s (now n)
         POO = SQRT ( ( ELAB + AM(KPROJ) ) * ( ELAB - AM(KPROJ) ) )
         PTTOT  = POO
         PXTTOT = PTTOT*TXI
         PYTTOT = PTTOT*TYI
         PZTTOT = PTTOT*TZI
      GO TO 2999
*  |
*  +-------------------------------------------------------------------*
*  |  AXsi0s are considered as nbars: after the interaction
*  |  a possible way to recover strangeness and charge conservation
*  |  is to flip an eventual nbar into a AXsi0,
*  |  or a pbar --> ASigma- and a pi+ --> K+ (or a pi0 --> K0) etc.
*  |  (we must flip two dbar quarks into sbar quarks)
 2035 CONTINUE
         LCHTYP = .TRUE.
         KPROJ = 9
*  |   Adjusting mass and momentum for AXsi0s (now nbar)
         POO = SQRT ( ( ELAB + AM(KPROJ) ) * ( ELAB - AM(KPROJ) ) )
         PTTOT  = POO
         PXTTOT = PTTOT*TXI
         PYTTOT = PTTOT*TYI
         PZTTOT = PTTOT*TZI
      GO TO 2999
*  |
*  +-------------------------------------------------------------------*
*  |  Xsi-s are considered as protons: after the interaction
*  |  a possible way to recover strangeness and charge conservation
*  |  is to flip an eventual proton into a Xsi-,
*  |  or a neutron --> Sigma- and a pi+ --> K0bar (or a pi0 --> K-) etc.
*  |  (we must flip two u quarks into s quarks)
 2036 CONTINUE
         LCHTYP = .TRUE.
         KPROJ = 1
*  |   The following card must be commented following the strangeness
*  |   and charge conservation recovery implementation
         ICHTOT = ICHTOT + 2
*  |   Adjusting mass and momentum for Xsi-s (now p)
         POO = SQRT ( ( ELAB + AM(KPROJ) ) * ( ELAB - AM(KPROJ) ) )
         PTTOT  = POO
         PXTTOT = PTTOT*TXI
         PYTTOT = PTTOT*TYI
         PZTTOT = PTTOT*TZI
      GO TO 2999
*  |
*  +-------------------------------------------------------------------*
*  |  AXsi+s are considered as pbars: after the interaction
*  |  a possible way to recover strangeness and charge conservation
*  |  is to flip an eventual pbar into a AXsi+,
*  |  or a nbar --> ASigma+ and a pi- --> K0 (or a pi0 --> K+) etc.
*  |  (we must flip two ubar quarks into sbar quarks)
 2037 CONTINUE
         LCHTYP = .TRUE.
         KPROJ = 2
*  |   The following card must be commented following the strangeness
*  |   and charge conservation recovery implementation
         ICHTOT = ICHTOT - 2
*  |   Adjusting mass and momentum for AXsi+s (now pbar)
         POO = SQRT ( ( ELAB + AM(KPROJ) ) * ( ELAB - AM(KPROJ) ) )
         PTTOT  = POO
         PXTTOT = PTTOT*TXI
         PYTTOT = PTTOT*TYI
         PZTTOT = PTTOT*TZI
      GO TO 2999
*  |
*  +-------------------------------------------------------------------*
*  |  Omega-s are considered as protons: after the interaction there
*  |  are many possible ways to recover strangeness and charge
*  |  conservation (we must flip a d quark into an s quark and two u
*  |  quarks into s quarks).
 2038 CONTINUE
         LCHTYP = .TRUE.
         KPROJ = 1
*  |   The following card must be commented following the strangeness
*  |   and charge conservation recovery implementation
         ICHTOT = ICHTOT + 2
*  |   Adjusting mass and momentum for Omega-s (now p)
         POO = SQRT ( ( ELAB + AM(KPROJ) ) * ( ELAB - AM(KPROJ) ) )
         PTTOT  = POO
         PXTTOT = PTTOT*TXI
         PYTTOT = PTTOT*TYI
         PZTTOT = PTTOT*TZI
      GO TO 2999
*  |
*  +-------------------------------------------------------------------*
*  |  Omegabar+s are considered as pbars: after the interaction there
*  |  are many possible ways to recover strangeness and charge
*  |  conservation (we must flip a dbar quark into an sbar quark and
*  |  two ubar quarks into sbar quarks).
 2039 CONTINUE
         LCHTYP = .TRUE.
         KPROJ = 2
*  |   The following card must be commented following the strangeness
*  |   and charge conservation recovery implementation
         ICHTOT = ICHTOT - 2
*  |   Adjusting mass and momentum for AOmega+s (now pbar)
         POO = SQRT ( ( ELAB + AM(KPROJ) ) * ( ELAB - AM(KPROJ) ) )
         PTTOT  = POO
         PXTTOT = PTTOT*TXI
         PYTTOT = PTTOT*TYI
         PZTTOT = PTTOT*TZI
      GO TO 2999
*  |
*  +-------------------------------------------------------------------*
 2999 CONTINUE
*  +-------------------------------------------------------------------*
*  |                                   Check for hydrogen!!!!
      IF ( IBTAR .EQ. 1 ) THEN
         TV = 0.D+00
         NP = NP0
         ITTA = 8 - 7 * ICHTAR
         IF ( KPROJ .EQ. 14 .AND. EKE .LT. 2.D+00 * GAMMIN ) THEN
            CALL PMPRAB ( KPROJ, EKE, POO, TXI, TYI, TZI, WE )
            RETURN
         END IF
*  |   Set a flag to avoid elastic collision in Hadrin
         IELFLG = +1
*  |   Set a flag to fully exploit charge exchange
         ICXFLG = 0
         CALL HADRIV ( KPROJ, POO, ELAB, TXI, TYI, TZI, ITTA )
         IELFLG = -1
         ICXFLG = -1
*  |  +----------------------------------------------------------------*
*  |  |
         IF ( IH .LE. 1 ) THEN
*  |  |  +-------------------------------------------------------------*
*  |  |  |
            IF ( ITH (1) .EQ. KPROJ ) THEN
               IH  = IH + 1
               ITH (2) = 1
               CXH (2) = TXI
               CYH (2) = TYI
               CZH (2) = TZI
               PLH (2) = 0.D+00
               ELH (2) = AM (1)
*  |  |  |
*  |  |  +-------------------------------------------------------------*
*  |  |  |
            ELSE IF ( ITH (1) .EQ. 1 ) THEN
               IH  = IH + 1
               ITH (2) = KPROJ
               CXH (2) = TXI
               CYH (2) = TYI
               CZH (2) = TZI
               PLH (2) = 0.D+00
               ELH (2) = AM (KPROJ)
*  |  |  |
*  |  |  +-------------------------------------------------------------*
*  |  |  |
            ELSE
            END IF
*  |  |  |
*  |  |  +-------------------------------------------------------------*
          END IF
*  |  |
*  |  +----------------------------------------------------------------*
 
*  |  +----------------------------------------------------------------*
*  |  |                  Looping over the particles produced in hadrin
         DO 110 J=1,IH
            NP = NP + 1
            TKI(NP) = ELH(J) - AM(ITH(J))
*  |  |  +-------------------------------------------------------------*
*  |  |  |
            IF ( TKI(NP) .GE. 0.D+00 ) THEN
               KPART(NP) = ITH(J)
               CXR(NP)   = CXH(J)
               CYR(NP)   = CYH(J)
               CZR(NP)   = CZH(J)
               PLR(NP)   = PLH(J)
               WEI(NP)   = WE
*  |  |  |  Updating conservation counters
               ENUCR  = ENUCR + ELH(J)
               PXNUCR = PXNUCR + PLR(NP)*CXR(NP)
               PYNUCR = PYNUCR + PLR(NP)*CYR(NP)
               PZNUCR = PZNUCR + PLR(NP)*CZR(NP)
               ICNUCR = ICNUCR + ICH(KPART(NP))
               IBNUCR = IBNUCR + IBAR(KPART(NP))
*  |  |  |
*  |  |  +-------------------------------------------------------------*
*  |  |  |
            ELSE IF ( ABS ( TKI (NP) ) .LE. ANGLGB ) THEN
*  |  |  |  Updating conservation counters
               ENUCR  = ENUCR  + ELH(J)
               ICNUCR = ICNUCR + ICH(ITH(J))
               IBNUCR = IBNUCR + IBAR(ITH(J))
               NP = NP - 1
*  |  |  |
*  |  |  +-------------------------------------------------------------*
*  |  |  |
            ELSE
               NP = NP - 1
            END IF
*  |  |  |
*  |  |  +-------------------------------------------------------------*
 110     CONTINUE
*  |  |
*  |  +----------------------------------------------------------------*
         KTARP = 1
         ANOW = 0.D+00
         ZNOW = 0.D+00
         ICRES = 0
         IBRES = 0
         ERES  = 0.D+00
         EOTEST = EOTEST - ENUCR
*  |  +----------------------------------------------------------------*
*  |  |
         IF ( ABS (EOTEST) .GT. ETEPS ) THEN
            WRITE (LUNERR,*)' Eventq: eotest failure with Hadrin',
     &                        EOTEST
            LRESMP = .TRUE.
            GO TO 50
*  |
*  +-<|--<--<--<--<--< go to resampling
         END IF
*  |  |
*  |  +----------------------------------------------------------------*
         RETURN
      END IF
*  |                                       end hydrogen check
*  +-------------------------------------------------------------------*
*
*  Now calling Nucrin to produce secondary particles, which are put
*  into the /Finuc/ common
*
      CALL NUCRIV ( KPROJ, ELAB, TXI, TYI, TZI, BBTAR,
     &              ZZTAR, RHO(IMAT) )
      IF ( LRESMP ) GO TO 50
*--<--<--<--<--< go to resampling if something was wrong
      TVGRE0 = TV - TVGREY
      TVEUZ  = 0.D+00
      TVTENT = TV
      TV = 0.D+00
      IF (NP .GT. MXP) GO TO 400
*  +-------------------------------------------------------------------*
*  |                     Looping over the particles produced in nucrin
*  |                     No need for Nucrin produced particles to
*  |                     change the numbering
      DO 101 I=NP0+1,NP
         I1 = KPART(I)
         WEI(I) = WE
         TKI(I) = TKI(I) - AM(I1)
         TXYZ   = CXR (I)**2 + CYR (I)**2 + CZR (I)**2
         IF ( ABS (TXYZ-1.D0) .GT. CSNPRN ) THEN
            WRITE ( LUNERR,* )
     &      ' **** Bad normalized cosines from Nucriv!!!',I,TXYZ,CXR(I),
     &        CYR(I),CZR(I)
            TXYZ = TXYZ - 1.D+00
            TXYZ = 1.D+00 - 0.5D+00 * TXYZ + 0.375D+00 * TXYZ * TXYZ
            CXR (I) = CXR (I) * TXYZ
            CYR (I) = CYR (I) * TXYZ
            CZR (I) = CZR (I) * TXYZ
            TXYZ   = CXR (I)**2 + CYR (I)**2 + CZR (I)**2
            IF ( ABS (TXYZ-1.D0) .GT. CSNPRN ) THEN
               LRESMP = .TRUE.
               GO TO 50
            END IF
         END IF
         IF (TKI(I) .GE. 0.D+00) GO TO 101
            TKI(I) = 0.D+00
            PLR(I) = 0.D+00
 101  CONTINUE
*  |
*  +-------------------------------------------------------------------*
      GO TO 500
C
C********************************************************************
C     GENERATE FIRST THE LOW ENERGY NUCLEONS FROM THE INTRANUCLEAR
C     CASCADE - CHECK IF ACCEPTABLE.
C     NOTE THAT AINEL AND RAKEKA ARE USING THE OLD PARTICLE NUMBERING
C********************************************************************
C
 200  CONTINUE
      NP = NP0
      LNUCRI = .FALSE.
*  +-------------------------------------------------------------------*
*  |                                   Check for hydrogen!!!!
      IF ( IBTAR .EQ. 1 ) THEN
         ANCOLL = 1.D+00
         NHAD = 0
         KPROJ = IJ
         KTARG = 1
         EPROJ = EKE + AM (KPROJ)
         UMO   = SQRT ( AM (1)**2 + AM (KPROJ)**2 + 2.D+00 * AM (1) *
     &                  EPROJ )
*  |  +----------------------------------------------------------------*
*  |  | Now diffractive events are taken into account for projectile
*  |  | -proton collisions!!! A. Ferrari and J. Ranft, 25-3-90
         IF ( LDIFFR (KPTOIP(KPROJ)) ) THEN
            CALL GRNDM(RNDM,1)
            IF ( RNDM (1) .LT. FRDIFF ) THEN
               IF ( POO .GT. PTHDFF ) THEN
                  LEVDIF = .TRUE.
                  CALL DIFEVV ( NHAD, KPROJ, KTARG, POO, EPROJ, UMO )
               ELSE
                  CALL HADEVV ( NHAD, KPROJ, KTARG, POO, EPROJ, UMO )
               END IF
            ELSE
               CALL HADEVV ( NHAD, KPROJ, KTARG, POO, EPROJ, UMO )
            END IF
*  |  |
*  |  +----------------------------------------------------------------*
*  |  |
         ELSE
            CALL HADEVV ( NHAD, KPROJ, KTARG, POO, EPROJ, UMO )
         END IF
*  |  |
*  |  +----------------------------------------------------------------*
*  |  +----------------------------------------------------------------*
*  |  |                   Loop over particles produced in hadevt/difevt
         DO 207 J=1,NHAD
            NP  = NP + 1
            CFE = 1.D+00
            SFE = 0.D+00
*  |  |   Get the "paprop" index
            KPART(NP) = KPTOIP (NREH(J))
            IF ( KPART (NP) .EQ. 0 ) THEN
               WRITE (LUNERR,*)
     &         ' **** Charmed particle produced in Hadevv/Difevv',
     &         NREH(J),HEPH(J),AMH(J),' ****'
               KPART (NP) = NREH (J)
            END IF
            PLR(NP) = SQRT ( PXH(J)**2 + PYH(J)**2 + PZH(J)**2 )
            PTH = SQRT( PXH(J)**2 + PYH(J)**2 )
            CDE = PZH(J) / PLR(NP)
            SDE = PTH    / PLR(NP)
            IF (SDE .GE. ANGLGB) THEN
               CFE = PXH(J) / PTH
               SFE = PYH(J) / PTH
            END IF
            CALL TTRANS ( TXI, TYI, TZI, CDE, SDE, SFE, CFE,
     &                    CXR(NP), CYR(NP), CZR(NP) )
            TKI(NP) = HEPH(J) - AMH(J)
            WEI(NP) = WE
*  |  |  +-------------------------------------------------------------*
*  |  |  |                   Updating conservation counters
            IF ( TKI(NP) .GT. 0.D+00 ) THEN
               EUZ = EUZ + HEPH(J)
               PUX = PUX + PLR(NP)*CXR(NP)
               PUY = PUY + PLR(NP)*CYR(NP)
               PUZ = PUZ + PLR(NP)*CZR(NP)
               ICU = ICU + ICH(NREH(J))
               IBU = IBU + IBAR(NREH(J))
*  |  |  |
*  |  |  +-------------------------------------------------------------*
*  |  |  |
            ELSE
               NP = NP-1
            END IF
*  |  |  |
*  |  |  +-------------------------------------------------------------*
 207     CONTINUE
*  |  |                                end loop
*  |  +----------------------------------------------------------------*
         ANOW = 0.D+00
         ZNOW = 0.D+00
         KTARP = 1
         ICRES = 0
         IBRES = 0
         ERES  = 0.D+00
         EOTEST = EOTEST - EUZ
*  |  +----------------------------------------------------------------*
*  |  |
         IF ( ABS (EOTEST) .GT. ETEPS ) THEN
            WRITE (LUNERR,*)' Eventq: eotest failure with',
     &                      ' Hadevt/Diffevt',EOTEST
         END IF
*  |  |
*  |  +----------------------------------------------------------------*
         RETURN
      END IF
*  |                                       end hydrogen check
*  +-------------------------------------------------------------------*
*  Now at first we sample the number of interactions
      CALL COREVT ( ZZTAR, BBTAR, IJ, POO, EKE )
*
*  First decide how much energy will go into the intranuclear cascade
*
      EKIN = EKE
*  Here starts the cascade particle production module!!!
      SOPPP (1) = EINCP
      SOPPP (2) = EINCN
      ARECL  = MAX ( ANOW - 1.D+00, 1.D+00 )
      KREJE  = 0
      NGREYT = NGREYP + NGREYN
      TVTENT = ANCOLL * ( AV0WEL - AEFRMA ) + TVGRE0
      TMP015 = 0.15D+00
      EEXTRA = 2.D+00* MIN ( TWOTWO * ( EKUPNU (1) + EKUPNU (2) )
     &                     , TMP015 *
     &       ( EKE - EINCP - EINCN - TVTENT ), EINCP + EINCN )
      EEXTRA = MAX ( EEXTRA, EKUPNU (1), EKUPNU (2) )
      LINCCK = .FALSE.
      PXCASC = 0.D+00
      PYCASC = 0.D+00
      PZCASC = 0.D+00
      PTXINT = 0.D+00
      PTYINT = 0.D+00
      PPINTR = 0.D+00
      EKRECL = 0.D+00
      IEXTRN = 0
      IEXTRP = 0
      LACCEP = .FALSE.
      IGREYT = 0
      IF ( NGREYT .EQ. 0 ) GO TO 206
*  +-------------------------------------------------------------------*
*  | Now looping until we reach the correct number of cascade
*  | nucleons
 205  CONTINUE
         IGREYT = IGREYP + IGREYN
         DSOPP  = SOPPP (1) + SOPPP (2)
*  |  Now it is not possible to reduce too heavily the available energy
         IF ( EKIN - TVGRE0 .LT. 0.5D+00 * EKE .OR.  PZCASC .GE.
     &        0.5D+00 * POO ) GO TO 206
*  |  +----------------------------------------------------------------*
*  |  |
         IF ( IGREYT .GE. NGREYT ) THEN
            AEXTRA = 2.D+00 * DSOPP / ( EKINAV (1) + EKINAV (2)
     &             + ESWELL (1) + ESWELL (2) )
     &             / ( IGREYT - NGREYT + 1 )
            CALL GRNDM(RNDM,1)
            RNDUMO = RNDM (1)
*  |  |  +-------------------------------------------------------------*
*  |  |  |
            IF ( RNDUMO .LT. AEXTRA .AND. LACCEP .AND. NINT (ANOW)
     &           .GT. 0 ) THEN
               REJE = ZNOW / ANOW
*  |  |  |  +----------------------------------------------------------*
*  |  |  |  |
               CALL GRNDM(RNDM,1)
               IF ( RNDM (1) .LT. REJE ) THEN
                  IEXTRP = MIN ( 1, NINT (ZNOW), NGREYP / 2 )
                  IF ( NGREYP + IEXTRP - IGREYP .LE. 0 ) GO TO 206
*  |  |  |  |
*  |  |  |  +----------------------------------------------------------*
*  |  |  |  |
               ELSE
*                 IEXTRN = MIN ( IEXTRN + 1, NINT (ANOW - ZNOW), 2 )
                  IEXTRN = MIN ( 1, NINT (ANOW - ZNOW), NGREYN / 2 )
                  IF ( NGREYN + IEXTRN - IGREYN .LE. 0 ) GO TO 206
               END IF
*  |  |  |  |
*  |  |  |  +----------------------------------------------------------*
*  |  |  |
*  |  |  +-------------------------------------------------------------*
*  |  |  |
            ELSE
               GO TO 206
*  |  |  |
*  |  |  |-->-->-->-->--> We have finished particles!!
            END IF
*  |  |  |
*  |  |  +-------------------------------------------------------------*
*  |  |
*  |  +----------------------------------------------------------------*
*  |  |
         ELSE IF ( DSOPP .LE. - EEXTRA ) THEN
            GO TO 206
*  |  |
*  |  |-->-->-->-->--> We have finished energy
         END IF
*  |  |
*  |  +----------------------------------------------------------------*
         LACCEP = .FALSE.
*  |  +----------------------------------------------------------------*
*  |  | Check how many nucleons are still available for cascade
         IF ( ANOW .LT. 1.5D+00 ) THEN
*  |  |  +-------------------------------------------------------------*
*  |  |  | Check if at least two high energy collisions are foreseen
*  |  |  | if yes the proper energy / momentum balance for the last
*  |  |  | two nucleons will be done inside Ferevv
            IF ( ANCOLL .GT. 1.5D+00 ) THEN
*  |  |  |
*  |  |  +-------------------------------------------------------------*
*  |  |  |     Only the valence collision is foreseen: we are forced
*  |  |  |     to perform here the balance for the 2 last nucleons!
*  |  |  |     (One might be used here for cascading, the other
*  |  |  |     will be used by Hadevv / Difevv)
            ELSE
               LLASTN = .TRUE.
               KTLAST = IJTARG (1)
               AMLAST = AM ( KTLAST )
*  |  |  |  +----------------------------------------------------------*
*  |  |  |  |
               IF ( NINT ( ZNOW ) .GT. 0 ) THEN
                  KTINC = 1
                  IGREYP = IGREYP + 1
*  |  |  |  |
*  |  |  |  +----------------------------------------------------------*
*  |  |  |  |
               ELSE
                  KTINC = 8
                  IGREYN = IGREYN + 1
               END IF
*  |  |  |  |
*  |  |  |  +----------------------------------------------------------*
               AMINC = AM (KTINC)
               UMIN2 = ( AMLAST + AMINC )**2
               DSOPP = MAX ( DSOPP, AV0WEL - AEFRMA )
  216          CONTINUE
               EKIN   = EKIN  - DSOPP
*  |  |  |  +----------------------------------------------------------*
*  |  |  |  |  Check if we can get enough invariant mass
               IF ( EKIN .LE. 0.2D+00 * EKE ) THEN
                  WRITE ( LUNOUT, * )
     &            ' *** Eventv: impossible to get enough invariant',
     &            ' mass for the last two nucleons! ****'
                  WRITE ( LUNERR, * )
     &            ' *** Eventv: impossible to get enough invariant',
     &            ' mass for the last two nucleons! ****'
                  LRESMP = .TRUE.
                  GO TO 50
*  |  |  |  |---<---<---<  Go to resampling
               END IF
*  |  |  |  |
*  |  |  |  +----------------------------------------------------------*
               ELEFT  = ETTOT - EINTR  - EKIN - AM (IJ)
               PLA    = SQRT ( EKIN * ( EKIN + 2.D+00 * AM (IJ) ) )
               PXLEFT = PXTTOT - PXINTR - PLA * TXI
               PYLEFT = PYTTOT - PYINTR - PLA * TYI
               PZLEFT = PZTTOT - PZINTR - PLA * TZI
*  |  |  |  Now we will divide Eleft and Pleft between the two
*  |  |  |  left nucleons!
               UMO2 = ELEFT**2 - PXLEFT**2 - PYLEFT**2 - PZLEFT**2
*  |  |  |  +----------------------------------------------------------*
*  |  |  |  |
               IF ( UMO2 .LE. UMIN2 ) THEN
*  |  |  |  |  +-------------------------------------------------------*
*  |  |  |  |  |
                  IF ( KTINC .EQ. 1 ) THEN
                     EINCP = EINCP + DSOPP
*  |  |  |  |  |
*  |  |  |  |  +-------------------------------------------------------*
*  |  |  |  |  |
                  ELSE
                     EINCN = EINCN + DSOPP
                  END IF
*  |  |  |  |  |
*  |  |  |  |  +-------------------------------------------------------*
                  GO TO 216
               END IF
*  |  |  |  |
*  |  |  |  +----------------------------------------------------------*
               UMO = SQRT ( UMO2 )
               ELACMS = 0.5D+00 * ( UMO2 + AMLAST**2 - AMINC**2 )
     &                / UMO
               EINCMS = UMO - ELACMS
               PCMS   = SQRT ( ELACMS**2 - AMLAST**2 )
               GAMCM = ELEFT  / UMO
               ETAX  = PXLEFT / UMO
               ETAY  = PYLEFT / UMO
               ETAZ  = PZLEFT / UMO
               CALL RACO ( CXXINC, CYYINC, CZZINC )
               PCMSX = PCMS * CXXINC
               PCMSY = PCMS * CYYINC
               PCMSZ = PCMS * CZZINC
*  |  |  |  Now go back from the CMS frame to the lab frame!!!
*  |  |  |  First the inc nucleon:
               ETAPCM = PCMSX * ETAX + PCMSY * ETAY + PCMSZ * ETAZ
               NP = NP + 1
               IF (NP .GT. MXP) GO TO 400
               KPART(NP)= KTINC
               WEI (NP) = WE
               TKI (NP) = GAMCM * EINCMS + ETAPCM - AMINC
               PHELP  = ETAPCM / (GAMCM + 1.D0) + EINCMS
               CXXINC = PCMSX + ETAX * PHELP
               CYYINC = PCMSY + ETAY * PHELP
               CZZINC = PCMSZ + ETAZ * PHELP
*  |  |  |  Updating conservation counters
               EINTR  = EINTR  + TKI(NP) + AMINC
               PXINTR = PXINTR + CXXINC
               PYINTR = PYINTR + CYYINC
               PZINTR = PZINTR + CZZINC
               ICINTR = ICINTR + ICH  (KPART(NP))
               IBINTR = IBINTR + IBAR (KPART(NP))
               PLR (NP) = SQRT (CXXINC * CXXINC + CYYINC * CYYINC +
     &                          CZZINC * CZZINC)
               CXR (NP) = CXXINC / PLR (NP)
               CYR (NP) = CYYINC / PLR (NP)
               CZR (NP) = CZZINC / PLR (NP)
*  |  |  |  Now the high energy collision nucleon
               EKLAST = GAMCM * ELACMS - ETAPCM - AMLAST
               PHELP  = - ETAPCM / (GAMCM + 1.D0) + ELACMS
               PXLAST = - PCMSX + ETAX * PHELP
               PYLAST = - PCMSY + ETAY * PHELP
               PZLAST = - PCMSZ + ETAZ * PHELP
*  |  |  |  Now we perform a Lorentz transformation to the rest fra-
*  |  |  |  me of the "last" nucleon, with the projectile momentum
*  |  |  |  along the +z axis
               AMPROJ = AM (IJ)
               EPROJ  = EKIN + AMPROJ
               ECHCK  = EPROJ + AMLAST + EKLAST
               PXCHCK = PLA * TXI + PXLAST
               PYCHCK = PLA * TYI + PYLAST
               PZCHCK = PLA * TZI + PZLAST
               UMO  = SQRT ( ECHCK**2 - PXCHCK**2 - PYCHCK**2 -
     &                       PZCHCK**2 )
               EPROJX = 0.5D+00 * ( UMO**2 - AMPROJ**2 - AMLAST**2 )
     &                / AMLAST
               PPROJX = SQRT ( EPROJX**2 - AMPROJ**2 )
               ETOTX  = EPROJX + AMLAST
*  |  |  |  Now set the parameters for the Lorentz transformation
               AAFACT = ECHCK  + ETOTX
               BBFACT = PPROJX - PZCHCK
               DDENOM = ETOTX * AAFACT - PPROJX * BBFACT
               GAMTRA = ( ECHCK * AAFACT + PPROJX * BBFACT ) / DDENOM
               ETAZ = - BBFACT * AAFACT / DDENOM
               ETAX = PXCHCK * ( GAMTRA + 1.D+00 ) / AAFACT
               ETAY = PYCHCK * ( GAMTRA + 1.D+00 ) / AAFACT
               PLABS = PPROJX
               ELABS = EPROJX
*  |  |  |  +----------------------------------------------------------*
*  |  |  |  |
               IF ( .NOT. LDIFFR (KPTOIP(IJ)) .OR. PLABS .LE. PTHDFF )
     &            THEN
                  RUUN = 1.D0
*  |  |  |  |
*  |  |  |  +----------------------------------------------------------*
*  |  |  |  |
               ELSE
                  CALL GRNDM(RNDM,1)
                  RUUN = RNDM (1)
               END IF
*  |  |  |  |
*  |  |  |  +----------------------------------------------------------*
               NHAD = 0
*  |  |  |  +----------------------------------------------------------*
*  |  |  |  |   Diffractive event
               IF ( RUUN .LE. FRDIFF ) THEN
                  LEVDIF = .TRUE.
                  CALL DIFEVV ( NHAD, IJ, KTLAST, PLABS, ELABS, UMO )
*  |  |  |  |
*  |  |  |  +----------------------------------------------------------*
*  |  |  |  |   Usual event
               ELSE
                  CALL HADEVV ( NHAD, IJ, KTLAST, PLABS, ELABS, UMO )
               END IF
*  |  |  |  |
*  |  |  |  +----------------------------------------------------------*
*  |  |  |  Now we have the particles in the nucleon rest frame,
*  |  |  |  transform back in the lab frame
*  |  |  |  +----------------------------------------------------------*
*  |  |  |  |           Looping over the produced particles
               DO 236 I = 1, NHAD
                  NP = NP + 1
                  IF (NP .GT. MXP) GO TO 400
                  CALL ALTRA ( GAMTRA, ETAX, ETAY, ETAZ, PXH (I),
     &                         PYH (I), PZH (I), HEPH (I), PLR (NP),
     &                         CXR (NP), CYR (NP), CZR (NP), ELR )
*  |  |  |  |  Updating conservation counters
                  EUZ  = EUZ + ELR
                  PUX  = PUX + CXR (NP)
                  PUY  = PUY + CYR (NP)
                  PUZ  = PUZ + CZR (NP)
                  KPART (NP) = KPTOIP ( NREH (I) )
                  IF ( KPART (NP) .EQ. 0 ) THEN
                     WRITE (LUNERR,*)
     &               ' **** Charmed particle produced in Hadevv',
     &               NREH(I),ELR,AMH(I),' ****'
                     KPART (NP) = NREH (I)
                  END IF
                  ICU = ICU + ICH  (NREH(I))
                  IBU = IBU + IBAR (NREH(I))
                  CXR (NP) = CXR (NP) / PLR (NP)
                  CYR (NP) = CYR (NP) / PLR (NP)
                  CZR (NP) = CZR (NP) / PLR (NP)
                  TKI (NP) = ELR - AMH (I)
                  WEI (NP) = WE
  236          CONTINUE
*  |  |  |  |
*  |  |  |  +----------------------------------------------------------*
*  |  |  | Now perform the standard tests for energy and momentum
*  |  |  | conservation
               ERES  = ETTOT  - EUZ - ENUCR  - EINTR
               PXRES = PXTTOT - PUX - PXNUCR - PXINTR
               PYRES = PYTTOT - PUY - PYNUCR - PYINTR
               PZRES = PZTTOT - PUZ - PZNUCR - PZINTR
               ICRES = ICHTOT - ICU - ICNUCR - ICINTR
               IBRES = IBTOT  - IBU - IBNUCR - IBINTR
               PTRES = MAX ( ABS (PXRES), ABS (PYRES), ABS (PZRES) )
*  |  |  |  +----------------------------------------------------------*
*  |  |  |  |
               IF ( IBRES .NE. 0 .OR. ICRES .NE. 0 .OR. ABS (ERES)
     &              .GT. 1.D-10*EPROJ .OR. PTRES .GT. 1.D-10*PTTOT )
     &              THEN
                  WRITE ( LUNERR, * )
     &                   ' Eventq: last nucleon failure!!', ICRES,
     &                     IBRES, REAL (ERES), REAL (PXRES),
     &                     REAL (PYRES), REAL (PZRES)
                  LRESMP = .TRUE.
                  GO TO 50
*  |  |  |--|--<--<--<--<--< go to resampling
               END IF
*  |  |  |  |
*  |  |  |  +----------------------------------------------------------*
               PTRES = 0.D+00
               RETURN
            END IF
*  |  |  |
*  |  |  +-------------------------------------------------------------*
         END IF
*  |  |
*  |  +----------------------------------------------------------------*
 
*  |  +----------------------------------------------------------------*
*  |  |
         IF ( NINT ( ANOW - ZNOW ) .LE. 0 ) THEN
            REJE = 1.D+00
*  |  |
*  |  +----------------------------------------------------------------*
*  |  |
         ELSE
            REJE = MAX ( ZNOW * ( NGREYP + IEXTRP - IGREYP ),
     &                   ZERZER )
            HELP = REJE + ( ANOW - ZNOW ) * ( NGREYN +
     &             IEXTRN - IGREYN )
            IF ( HELP .GT. 0.D+00 ) THEN
               REJE = REJE / HELP
            ELSE
               REJE = 1.D+00
               WRITE(LUNERR,*)' EVENTV: HELP=0, ANOW,ZNOW',ANOW,ZNOW
     &                       ,' NGREYN,IGREYN,IEXTRN',
     &                          NGREYN,IGREYN,IEXTRN,
     &                        ' NGREYP,IGREYP,IEXTRP',
     &                          NGREYP,IGREYP,IEXTRP
            END IF
         END IF
*  |  |
*  |  +----------------------------------------------------------------*
 
*  |  +----------------------------------------------------------------*
*  |  |
         CALL GRNDM(RNDM,1)
         IF ( RNDM (1) .LT. REJE ) THEN
            ILO  = 1
            ILLO = 2
            KP   = 1
*  |  |
*  |  +----------------------------------------------------------------*
*  |  |
         ELSE
            ILO  = 2
            ILLO = 1
            KP   = 8
         END IF
*  |  |
*  |  +----------------------------------------------------------------*
         CALL RAKEKV ( ILO, EXSOP, BBTAR, TKIN, TSTRCK, PLA, ARECL,
     &                 TKRECL, EFERMI, CDE, SDE )
         IF ( EKIN - TKIN .LT. 0.5D+00 * EKE ) GO TO 206
*  |--|--<--<--<--<--< avoid to deplete too much the energy
*  |  +----------------------------------------------------------------*
*  |  |  Now check if the nucleon energy is acceptable:
         IF ( SOPPP (ILO) .LT. TKIN )  THEN
*  |  |  +-------------------------------------------------------------*
*  |  |  |
            IF ( TKIN - SOPPP (ILO) .LT. SOPPP (ILLO) + 1.5D+00 *
     &           EEXTRA ) THEN
               SOPPP (ILLO) = SOPPP (ILLO) + SOPPP (ILO) - TKIN
     &                      - EXSOP (ILO)
               SOPPP (ILO)  = 0.D+00
*  |  |  |
*  |  |  +-------------------------------------------------------------*
*  |  |  |
            ELSE
               KREJE = KREJE + 1
               IF ( KREJE .GT. 10 ) GO TO 206
               SOPPP (ILLO) = SOPPP (ILLO) + SOPPP (ILO)
               SOPPP (ILO)  = 0.D+00
               GO TO 205
            END IF
*  |  |  |
*  |  |  +-------------------------------------------------------------*
*  |  |
*  |  +----------------------------------------------------------------*
*  |  |
         ELSE
            SOPPP (ILO) = SOPPP (ILO) - TKIN - EXSOP (ILO)
         END IF
*  |  |
*  |  +----------------------------------------------------------------*
         IF ( NINT ( ANOW ) .LE. 0 ) GO TO 206
*  |  +----------------------------------------------------------------*
*  |  |   Update the current atomic and mass number
         IF ( ILO .EQ. 1 ) THEN
            IF ( NINT ( ZNOW ) .LE. 0 ) GO TO 206
            IGREYP = IGREYP + 1
            ZNOW   = ZNOW   - 1.D+00
*  |  |
*  |  +----------------------------------------------------------------*
*  |  |
         ELSE
            IGREYN = IGREYN + 1
         END IF
*  |  |
*  |  +----------------------------------------------------------------*
         ANOW = ANOW - 1.D+00
*  |  Now make the kinematical tests!!
         FRSURV = ARECL / BBTAR
         NP = NP + 1
         IF ( NP .GT. MXP ) GO TO 400
         EKIN   = EKIN - TKIN
*  Note the change!!
         LINCCK = ARECL .LT. 30.D+00 .AND. FRSURV .LE. 0.33D+00
*  |  +----------------------------------------------------------------*
*  |  |
         IF ( .NOT. LINCCK ) THEN
            CALL SFECFE ( SFE, CFE )
*  |  |
*  |  +----------------------------------------------------------------*
*  |  |
         ELSE
            IF ( CDE .GT. 0.D+00 ) THEN
               CDE = CDE * FRSURV / 0.33D+00
               SDE = SQRT ( 1.D+00 - CDE**2 )
            END IF
            CALL SFECFE ( SFE, CFE )
         END IF
*  |  |
*  |  +----------------------------------------------------------------*
         PTXINT = PTXINT + PLA * CFE * SDE
         PTYINT = PTYINT + PLA * SFE * SDE
         PPINTR = PPINTR + PLA * CDE
         CALL TTRANS ( TXI, TYI, TZI, CDE, SDE, SFE, CFE,
     &                 CXR (NP), CYR (NP), CZR (NP) )
*  | Generated nucleon is acceptable - save it
         LACCEP = .TRUE.
         TVGRE0 = TVGRE0 + EXSOP (ILO)
         TVGREY = TVGREY + TKIN + TKRECL - TSTRCK - EBNDNG (ILO)
         EKRECL = EKRECL + TKRECL
         ARECL  = MAX ( ARECL - 1.D+00, 1.D+00 )
         KPART(NP)= KP
         WEI(NP)  = WE
         TKI(NP)  = TSTRCK
         PLR(NP)  = PLA
*  |                        Updating the counters for spent momenta
*  |                        of the cascade nucleons (the momentum
*  |                        spent by the high energy particles)
         SINTH = SQRT ( 1.D+00 - COSTH * COSTH )
         PLA = SQRT ( ( EFERMI + TKIN )**2 - AMNUSQ (ILO) )
         PZLEFT = PLA * CDE - COSTH * PKFRMI
         IF ( PZLEFT + PZCASC .LT. EKE - EKIN + TVGRE0 ) THEN
            PZCASC = EKE - EKIN + TVGRE0
            PXCASC = PXCASC + CFE * ( PLA * SDE - SINTH * PKFRMI )
            PYCASC = PYCASC + SFE * ( PLA * SDE - SINTH * PKFRMI )
         ELSE
            PZCASC = PZCASC + PZLEFT
            PXCASC = PXCASC + CFE * ( PLA * SDE - SINTH * PKFRMI )
            PYCASC = PYCASC + SFE * ( PLA * SDE - SINTH * PKFRMI )
         END IF
*  |                        Updating conservation counters
         EINTR  = EINTR  + TKI(NP) + AM (KPART(NP))
         PXINTR = PXINTR + PLR(NP) * CXR(NP)
         PYINTR = PYINTR + PLR(NP) * CYR(NP)
         PZINTR = PZINTR + PLR(NP) * CZR(NP)
         ICINTR = ICINTR + ICH  (KPART(NP))
         IBINTR = IBINTR + IBAR (KPART(NP))
      GO TO 205
*  |
*  +--<--<--<--<--< go to sample another nucleon
 206  CONTINUE
*  +-------------------------------------------------------------------*
*  |  First check charge and baryonic number conservation:
      IF ( IGREYP .NE. ICINTR .OR. ( IGREYP + IGREYN ) .NE.
     &     IBINTR ) THEN
         WRITE ( LUNERR,* )' Eventq: charge or baryon number',
     &                     ' conservation failure in the Inc',
     &                     ' section!!', IGREYP, ICINTR,
     &                       IGREYP + IGREYN, IBINTR
         LRESMP = .TRUE.
         GO TO 50
*--|--<--<--<--<--< go to resampling
      END IF
*  |
*  +-------------------------------------------------------------------*
      TVTENT = TVGRE0
* Tvgrey is already inside Eincp and Eincn since Rakekv updates the
* estimated excitation energy (and Tvgre0, Eincp, Eincn and the
* Soppp array)
      EINCP  = EINCP - SOPPP (1)
      EINCN  = EINCN - SOPPP (2)
      EKIN   = EKIN  - TVTENT
*  +-------------------------------------------------------------------*
*  |   Check if we have not spent too much energy!!!
      IF ( EKIN .LT. 0.5D+00 * EKE ) THEN
         LRESMP = .TRUE.
         WRITE ( LUNERR,* ) ' Eventq: Ekin after inc too low!! ',
     &                        EKIN, EKE, IJIJ
      END IF
*  |
*  +-------------------------------------------------------------------*
      IF ( LRESMP ) GO TO 50
*
*--<--<--<--<--<--< go to resampling if something was wrong
*   Here the modification to take into account properly the
*   cascade nucleon momenta
*   New version: rotate the spent momentum in the lab frame
      IF ( IGREYT .GT. 0 ) THEN
         PZCASC = MAX ( PZCASC, ZERZER )
         PTH = PXCASC * PXCASC + PYCASC * PYCASC
         PLA = SQRT ( PTH + PZCASC * PZCASC )
         PTH = SQRT ( PTH )
         CDE = PZCASC / PLA
         SDE = PTH    / PLA
         IF (SDE .GE. ANGLGB) THEN
            CFE = PXCASC / PTH
            SFE = PYCASC / PTH
         ELSE
            CFE = 1.D+00
            SFE = 0.D+00
         END IF
         CALL TTRANS ( TXI, TYI, TZI, CDE, SDE, SFE, CFE,
     &                 CXXINC, CYYINC, CZZINC )
         PXCASC = PLA * CXXINC
         PYCASC = PLA * CYYINC
         PZCASC = PLA * CZZINC
      END IF
      PXLEFT = PXTTOT - PXCASC - TXI * TVGRE0
      PYLEFT = PYTTOT - PYCASC - TYI * TVGRE0
      PZLEFT = PZTTOT - PZCASC - TZI * TVGRE0
      PLA = SQRT ( PXLEFT * PXLEFT + PYLEFT * PYLEFT + PZLEFT
     &           * PZLEFT )
      DEKVSP = PLA - EKIN - AM (IJ)
*  +-------------------------------------------------------------------*
*  |  Check the momentum versus the total energy
      IF ( DEKVSP .GE. 0.D+00 ) THEN
*  |  +----------------------------------------------------------------*
*  |  |  There are good reasons to believe that the following attempt
*  |  |  is dangerous rather than useful, leading to a Dp > DEk
         IF ( TVGRE0 .GT. 0.D+00 ) THEN
            TVTENT = TVTENT - TVGRE0
            EKIN   = EKIN + TVGRE0
            TVGRE0 = 0.D+00
            PXLEFT = PXTTOT - PXCASC
            PYLEFT = PYTTOT - PYCASC
            PZLEFT = PZTTOT - PZCASC
            PLA = SQRT ( PXLEFT * PXLEFT + PYLEFT * PYLEFT + PZLEFT
     &                 * PZLEFT )
            IF ( EKIN + AM (IJ) .GT. PLA ) GO TO 300
*  |  |  This part is new (1-10-91)
            DEEXTR = EKE - EKIN - EINTR + IGREYP * AM (1)
     &             + IGREYN * AM (8) - 1.5D+00
     &             * (IGREYP+IGREYN) * EBNDAV
            IF ( DEEXTR .GT. 0.D+00 ) THEN
               EKIN  = EKIN  + 0.5D+00 * DEEXTR
               EINCP = EINCP - 0.5D+00 * DEEXTR
               EINCN = EINCN - 0.5D+00 * DEEXTR
               PXLEFT = PXTTOT - PXCASC + 0.5D+00 * TXI * DEEXTR
               PYLEFT = PYTTOT - PYCASC + 0.5D+00 * TYI * DEEXTR
               PZLEFT = PZTTOT - PZCASC + 0.5D+00 * TZI * DEEXTR
               PLA = SQRT ( PXLEFT * PXLEFT + PYLEFT * PYLEFT + PZLEFT
     &                    * PZLEFT )
               IF ( EKIN + AM (IJ) .GT. PLA ) GO TO 300
            END IF
         END IF
*  |  |
*  |  +----------------------------------------------------------------*
         IF ( PLA - EKIN - AM (IJ) .LE. 1.D-04 * PLA ) GO TO 300
*  |  +----------------------------------------------------------------*
*  |  |  Printing condition relaxed, A.F., 23-12-92
         IF ( PLA - EKIN - AM (IJ) .GT. 1.D-02 * PLA ) THEN
            WRITE ( LUNERR,* )' Eventv: ekin+am < pla,ij,igreyt',
     &                          EKIN+AM(IJ),PLA,IJ,IGREYT
         END IF
*  |  |
*  |  +----------------------------------------------------------------*
         PTH = PLA
         PLA = EKIN + AM(IJ)
         PXLEFT = PXLEFT * PLA / PTH
         PYLEFT = PYLEFT * PLA / PTH
         PZLEFT = PZLEFT * PLA / PTH
      END IF
*  |
*  +-------------------------------------------------------------------*
  300 CONTINUE
*  end new version
*  +-------------------------------------------------------------------*
*  |  Check if we have enough energy to enter the high energy module
*  |  if not reset the accumulators and go to nucrin
      IF ( PLA .LT. PTHRSH (IJIJ) ) THEN
         PREF = 0.45D+00 * POO
*  |  +----------------------------------------------------------------*
*  |  |  Check if the momentum is much smaller than the original one
         IF ( PLA .LE. PREF ) THEN
*  |  |  +-------------------------------------------------------------*
*  |  |  |
            IF ( PLA .LE. 0.4D+00 * POO ) THEN
               WRITE ( LUNERR,* )' Eventv: Pla < 0.4 Poo',PLA,POO,IJIJ,
     &                             TKIN,IMAT
               LRESMP = .TRUE.
               GO TO 50
            END IF
*  |  |  |
*  |  |  +-------------------------------------------------------------*
         END IF
*  |  |
*  |  +----------------------------------------------------------------*
         PREF = MAX ( PREF, TWOTWO )
         PREF = MIN ( PREF, FOUFOU )
*  |  +----------------------------------------------------------------*
*  |  |  Originally it was < 10, but with Fermi momentum sometimes
*  |  |  we were calling Hadrin with p > 10 GeV/c
*        IF ( POO .LE. 9.75D+00 ) THEN
         IF ( ( POO .LE. 9.75D+00 .AND. IJNUCR (IJIJ) .LE. 0 )
     &        .OR. PLA .LT. PREF ) THEN
            IGREYP = 0
            IGREYN = 0
            KTARP  = 0
            KTARN  = 0
            TVGRE0 = 0.D+00
            TVGREY = 0.D+00
            EINTR  = 0.D+00
            PXINTR = 0.D+00
            PYINTR = 0.D+00
            PZINTR = 0.D+00
            ICINTR = 0
            IBINTR = 0
            ANOW   = BBTAR
            ZNOW   = ZZTAR
            IF ( IJ .EQ. 26 ) IJ = 23
            GO TO 100
         END IF
*  |  |
*  |  +----------------------------------------------------------------*
      END IF
*  |
*  +-------------------------------------------------------------------*
      TXX = PXLEFT / PLA
      TYY = PYLEFT / PLA
      TZZ = PZLEFT / PLA
      ZTEMP = ZZTAR - IGREYP
      ATEMP = BBTAR - IGREYN
      ERES  = ETTOT - EKIN - AM (IJ) - EINTR
      AMNRES = AMUAMU * ATEMP + 1.D-03 * FKENER ( ATEMP, ZTEMP ) -
     &         ZTEMP * AMELEC + ELBNDE ( NINT (ZTEMP) )
C
C********************************************************************
C     FOR MOMENTA ABOVE 5.0 GEV/C USE NUCEVT
C********************************************************************
C
*
*  From here the high energy model....
*
      NNHAD = 0
      CALL NUCEVV ( NNHAD, IJ, PLA, EKIN, TXX, TYY, TZZ )
      IF ( LRESMP ) GO TO 50
*--<--<--<--<--< go to resampling if something was wrong
*  +-------------------------------------------------------------------*
*  |
      IF ( LLASTN .AND. NINT ( ANOW ) .EQ. 1 ) THEN
*  |  +----------------------------------------------------------------*
*  |  |
         IF ( KTINC .EQ. 1 ) THEN
            IGREYP = IGREYP + 1
            EINCP  = EINCP  + EKINC
            ZNOW   = ZNOW - 1.D+00
*  |  |
*  |  +----------------------------------------------------------------*
*  |  |
         ELSE
            IGREYN = IGREYN + 1
            EINCN  = EINCN  + EKINC
         END IF
*  |  |
*  |  +----------------------------------------------------------------*
         ANOW = ANOW - 1.D+00
         NP = NP + 1
         IF (NP .GT. MXP) GO TO 400
         KPART(NP)= KPTOIP (KTINC)
         WEI (NP) = WE
         TKI (NP) = EKINC
*  |  Updating conservation counters
         EINTR  = EINTR  + TKI (NP) + AMINC
         PXINTR = PXINTR + PXXINC
         PYINTR = PYINTR + PYYINC
         PZINTR = PZINTR + PZZINC
         ICINTR = ICINTR + ICH  (KTINC)
         IBINTR = IBINTR + IBAR (KTINC)
         PLR (NP) = SQRT (PXXINC * PXXINC + PYYINC * PYYINC +
     &                    PZZINC * PZZINC)
         CXR (NP) = PXXINC / PLR (NP)
         CYR (NP) = PYYINC / PLR (NP)
         CZR (NP) = PZZINC / PLR (NP)
      END IF
*  |
*  +-------------------------------------------------------------------*
 
*  +-------------------------------------------------------------------*
*  |                          Loop over particles produced in nucevt
      DO 301 I=1,NNHAD
         NP = NP+1
         IF (NP .GT. MXP) GO TO 400
*  |  Get the "paprop" index for Nucevt produced particles
         KPART(NP) = KPTOIP (NRENU(I))
         IF ( KPART (NP) .EQ. 0 ) THEN
            WRITE (LUNERR,*)' **** Charmed particle produced in Nucevv',
     &      NRENU(I),HEPNU(I),AMNU(I),' ****'
            KPART (NP) = NRENU (I)
         END IF
         PLR  (NP) = SQRT ( PXNU (I)**2 + PYNU (I)**2 + PZNU (I)**2 )
         CXR  (NP) = PXNU (I) / PLR (NP)
         CYR  (NP) = PYNU (I) / PLR (NP)
         CZR  (NP) = PZNU (I) / PLR (NP)
         TKI(NP) = HEPNU(I) - AMNU(I)
         WEI(NP) = WE
*  |  +----------------------------------------------------------------*
*  |  |
         IF (TKI(NP) .LE. 0.D+00) THEN
*  |  | Is this the only check on the generated particle energy??
            EUZ = EUZ - HEPNU(I)
            PUX = PUX - PXNU(I)
            PUY = PUY - PYNU(I)
            PUZ = PUZ - PZNU(I)
            ICU = ICU - ICHNU(I)
            IBU = IBU - IBARNU(I)
            WRITE (LUNERR,*) ' Eventq: Kin en. < 0 from nucevt',TKI(NP)
            NP  = NP - 1
         END IF
*  |  |
*  |  +----------------------------------------------------------------*
 301  CONTINUE
*  |
*  +-------------------------------------------------------------------*
 500  CONTINUE
*
*  Now we have to compute the excitation energy: we have used Eincp and
*  Eincn for cascade protons and neutrons, Tvgre0 and Tvgrey have been
*  generated by cascade nucleons under threshold, Tveuz has been gene-
*  rated by the high energy collisions: however Tveuz and Tvgrey are
*  only approximate since they have been computed starting from an
*  average binding energy, furthermore Tvgrey is already accounted for
*  inside Eincp and Eincn. The actual energy spent inside high energy
*  collisions was Eke - Eincp - Eincn - Tvgre0 + Efrm, the one in ca-
*  scade nucleons (approximately) Eincp + Eincn - Tvgrey: so first
*  check the energy balance without excitation energy
*  and then compute Tv!!!
*
*  Now the balance!!!!!!!!!!!!!!!
*
      ERES  = ETTOT  - EUZ - ENUCR  - EINTR
      PXRES = PXTTOT - PUX - PXNUCR - PXINTR
      PYRES = PYTTOT - PUY - PYNUCR - PYINTR
      PZRES = PZTTOT - PUZ - PZNUCR - PZINTR
      ICRES = ICHTOT - ICU - ICNUCR - ICINTR
      IBRES = IBTOT  - IBU - IBNUCR - IBINTR
*  +-------------------------------------------------------------------*
*  |
      IF ( NINT (ZNOW) .NE. ICRES .OR. NINT (ANOW) .NE. IBRES ) THEN
*  |  +----------------------------------------------------------------*
*  |  |
         IF ( LNUCRI ) THEN
            WRITE (LUNERR,*)' Eventq: charge/baryon conservation',
     &                      ' failure with Nucrin',
     &                      NINT (ZNOW), ICRES, NINT (ANOW), IBRES
*  |  |
*  |  +----------------------------------------------------------------*
*  |  |
         ELSE
            WRITE (LUNERR,*)' Eventq: charge/baryon conservation',
     &                      ' failure with Nucevt',
     &                      NINT (ZNOW), ICRES, NINT (ANOW), IBRES
         END IF
*  |  |
*  |  +----------------------------------------------------------------*
         LRESMP = .TRUE.
         GO TO 50
*  |--<--<--<--<--< go to resampling
      END IF
*  |
*  +-------------------------------------------------------------------*
*  +-------------------------------------------------------------------*
*  |
      IF ( IBRES .GT. 0 ) THEN
         AMMRES = ANOW * AMUAMU + 0.001D+00 * FKENER ( ANOW, ZNOW )
         AMNRES = AMMRES - ZNOW * AMELEC + ELBNDE ( ICRES )
*        AMNRES = ANOW * AMUC12
         EKR0   = ERES - AMNRES
*  |  Now switch to atomic masses:
         ERES   = ERES + AMMTAR - AMNTAR - ( ZZTAR - ZNOW ) * AMELEC
         EKRES  = ERES - AMMRES
         TVTENT = TVGRE0 + TVGREY + TVEUZ
         TVCOMP = ERES - ANOW * AMUC12
         IF ( LNUCRI ) TVCOMP = TVCOMP + ( AMNTAR - BBTAR * AMUC12 )
*  |
*  +-------------------------------------------------------------------*
*  |
      ELSE
         AMMRES = 0.D+00
         AMNRES = 0.D+00
         ERES   = 0.D+00
         EKR0   = 0.D+00
         EKRES  = 0.D+00
         TVTENT = 0.D+00
         GO TO 600
      END IF
*  |
*  +-------------------------------------------------------------------*
*  Check that the kinetic energy of the residual nuclei is not much
*  different from our prevision and kinematically consistent with
*  its momentum
*  +-------------------------------------------------------------------*
*  |
      IF ( EKRES .LE. 0.D+00 ) THEN
         WRITE ( LUNERR,* )' Eventq: negative kinetic energy for',
     &                     ' the residual nucleus!!',ICRES,IBRES,
     &                       REAL (EKRES), LNUCRI
*  |  +----------------------------------------------------------------*
*  |  |
         IF ( EKRES .LT. -3.D-3 ) THEN
            LRESMP = .TRUE.
            GO TO 50
*--|--|--<--<--<--<--< go to resampling
         END IF
*  |  |
*  |  +----------------------------------------------------------------*
         EKRES  = 0.D+00
         TVRECL = 0.D+00
         AMSTAR = AMMRES
         TVCMS  = 0.D+00
         PTRES2 = 0.D+00
         PXRES  = 0.D+00
         PYRES  = 0.D+00
         PZRES  = 0.D+00
*  |
*  +-------------------------------------------------------------------*
*  |
      ELSE
         PTRES2 = PXRES**2 + PYRES**2 + PZRES**2
         AMSTAR = ERES**2  - PTRES2
*  |  +----------------------------------------------------------------*
*  |  |
         IF ( AMSTAR .GE. AMMRES**2 ) THEN
            AMSTAR = SQRT ( AMSTAR )
            TVCMS  = AMSTAR - AMMRES
*  |  |
*  |  +----------------------------------------------------------------*
*  |  |              If the following condition is satisfied it is only
*  |  |              a rounding problem, set the excitation energy to 0
*  |  |              and continue
         ELSE IF ( AMMRES**2 - AMSTAR .LT. 2.D+00 * AMSTAR * TVEPSI
     &             ) THEN
            AMSTAR = AMMRES
            TVCMS  = 0.D+00
*  |  |
*  |  +----------------------------------------------------------------*
*  |  |
         ELSE IF ( AMSTAR .LE. 0.D+00 ) THEN
            WRITE ( LUNERR,* )' Eventq: immaginary mass for',
     &                        ' the residual nucleus!!',ICRES,IBRES,
     &                          REAL (AMSTAR)
            LRESMP = .TRUE.
            GO TO 50
*--|--|--<--<--<--<--< go to resampling
*           AMSTAR = AMMRES
*           TVCMS  = 0.D+00
*  |  |
*  |  +----------------------------------------------------------------*
*  |  |
         ELSE
            AMSTAR = SQRT ( AMSTAR )
*  |  |  +-------------------------------------------------------------*
*  |  |  |           If the following condition is satisfied it is only
*  |  |  |           a rounding problem, set the excitation energy to 0
*  |  |  |           and continue
            IF ( AMMRES - AMSTAR .LT. TVEPSI ) THEN
               AMSTAR = AMMRES
               TVCMS  = 0.D+00
               TVRECL = ERES - AMSTAR
               GO TO 550
            END IF
*  |  |  |
*  |  |  +-------------------------------------------------------------*
            WRITE ( LUNERR,* )' Eventq: negative excitation energy for',
     &                        ' the residual nucleus!!',ICRES,IBRES,
     &                          REAL ( AMSTAR - AMMRES ), LNUCRI
*  |  |  +-------------------------------------------------------------*
*  |  |  |
            IF ( AMSTAR - AMMRES .LT. -3.D-3 ) THEN
               LRESMP = .TRUE.
               GO TO 50
*--|--|--|--<--<--<--<--< go to resampling
            END IF
*  |  |  |
*  |  |  +-------------------------------------------------------------*
            AMSTAR = AMMRES
            TVCMS  = 0.D+00
            AMSTAR = AMMRES
            TVCMS  = 0.D+00
            HELP   = SQRT ( ( ERES - AMMRES ) * ( ERES + AMMRES )
     &             / PTRES2 )
            PXRES = PXRES * HELP
            PYRES = PYRES * HELP
            PZRES = PZRES * HELP
            PTRES2 = PTRES2 * HELP * HELP
         END IF
*  |  |
*  |  +----------------------------------------------------------------*
         TVRECL = ERES - AMSTAR
      END IF
*  |
*  +-------------------------------------------------------------------*
  550 CONTINUE
*  The following two cards are equivalent providing the kinematical
*  limits are ok and we use for both amnres or ammres! Now Tv is left
*  = 0
      TV     = 0.D+00
*     TV     = EKRES
*     TV     = TVRECL + TVCMS
*  +-------------------------------------------------------------------*
*  |
      IF ( .NOT. LEVPRT .AND. ABS ( ( TVTENT - TVCOMP ) / TVCOMP )
     &     .GT. 30000000.D+00 ) THEN
*  |  +----------------------------------------------------------------*
*  |  |
         IF ( LINCTV ) THEN
            WRITE ( LUNERR,* )
     &                      ' Eventq: excitation energy very different',
     &                      ' from the approximate one!!', ICRES, IBRES,
     &                        REAL (TVTENT), REAL (TVCOMP), LNUCRI
         END IF
*  |  |
*  |  +----------------------------------------------------------------*
      END IF
*  |
*  +-------------------------------------------------------------------*
      EKRES = TVRECL
 600  CONTINUE
      EOTEST = EOTEST - EUZ - ENUCR - EINTR - EKR0 - AMNRES
*  +-------------------------------------------------------------------*
*  |
      IF ( ABS (EOTEST) .GT. ETEPS ) THEN
*  |  +----------------------------------------------------------------*
*  |  |
         IF ( LNUCRI ) THEN
            WRITE (LUNERR,*)' Eventq: eotest failure with Nucrin',
     &                      EOTEST
*  |  |
*  |  +----------------------------------------------------------------*
*  |  |
         ELSE
            WRITE (LUNERR,*)' Eventq: eotest failure with Nucevt',
     &                      EOTEST
         END IF
*  |  |
*  |  +----------------------------------------------------------------*
         LRESMP = .TRUE.
         GO TO 50
*  |--|--<--<--<--<--< go to resampling
      END IF
*  |
*  +-------------------------------------------------------------------*
      IF ( IBRES .EQ. 0 ) RETURN
*-->-->-->-->--> Here for the evaporation step!!!
*  +-------------------------------------------------------------------*
*  |                      Check if the evaporation is requested
      IF ( LEVPRT ) THEN
         PTRES = SQRT ( PTRES2 )
         CALL EVEVAP ( WE )
         IF ( LRESMP ) GO TO 50
*
*--|--<--<--<--<--<--< go to resampling
*  |
*  +-------------------------------------------------------------------*
*  |                      No evaporation
      ELSE
         TVHEAV = 0.D+00
      END IF
*  |
*  +-------------------------------------------------------------------*
      IF (IPRI.NE.1) RETURN
C
C********************************************************************
C     TEST PRINTOUT
C********************************************************************
C
      WRITE(LUNOUT,590)NP-NP0,NNHAD,IJ,IMAT,POO,EKE,TXI,TYI,TZI,WE
 590  FORMAT (4I7,6F12.6)
      DO 501 I=NP0+1,NP
      WRITE(LUNOUT,591)I,KPART(I),CXR(I),CYR(I),CZR (I),TKI(I),PLR (I),
     *                         WEI(I)
 591  FORMAT (2I5,6F12.6)
 501  CONTINUE
      RETURN
C
C********************************************************************
C     FINUC FLOWS OVER - THIS IS FATAL - INCREASE THE SIZE OF IT
C********************************************************************
C
 400  CONTINUE
      WRITE(LUNOUT,490)
 490  FORMAT(1X,'OVERFLOW IN EVENTQ - INCREASE THE SIZE OF THE',
     1' COMMON BLOCK FINUC.')
      STOP
      END
