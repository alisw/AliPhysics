*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:19:54  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.41  by  S.Giani
*-- Author :
*$ CREATE COREVT.FOR
*COPY COREVT
*                                                                      *
*=== corevt ===========================================================*
*                                                                      *
      SUBROUTINE COREVT ( ZZTAR, BBTAR, KPROJ, PPROJ, EKE )
 
#include "geant321/dblprc.inc"
#include "geant321/dimpar.inc"
#include "geant321/iounit.inc"
*
*----------------------------------------------------------------------*
*                                                                      *
*    Last change   on  07-may-92   by   Alfredo Ferrari, INFN-Milan    *
*                                                                      *
*    This version has been developed by A. Ferrari trying to take into *
*    account a few papers about correlations between projectile colli- *
*    sions and secondary collisions. An improved (but slower) version  *
*    is going on.                                                      *
*    This subroutine correlates intranuclear cascade with the number   *
*    of projectile collisions sampled in the nucleus                   *
*                                                                      *
*              input parameters:                                       *
*      bbtar - atomic number of the target                             *
*      kproj - type of the projectile                                  *
*      pproj - momentum of the projectile                              *
*                                                                      *
*              output parameters (in common corinc)                    *
*      frainc - reduction factor for intran. cascade energy,           *
*               inc. correlat.                                         *
*      nsea+1 - number of high energy collisions in the nucleus        *
*                                                                      *
*----------------------------------------------------------------------*
*
#include "geant321/balanc.inc"
#include "geant321/corinc.inc"
#include "geant321/nucdat.inc"
#include "geant321/parevt.inc"
#include "geant321/part.inc"
#include "geant321/qquark.inc"
#include "geant321/resnuc.inc"
*
      PARAMETER ( PIO2   = 0.5D+00 * PIPIPI )
      PARAMETER ( SQPI   = 1.772453850905516 D+00 )
      PARAMETER ( SQPIT3 = 3.D+00 * SQPI )
      PARAMETER ( TWOSQT = 2.828427124746190 D+00 )
      PARAMETER ( ECUTRF = 100.0 D+00 )
      PARAMETER ( UMOREF = 13.83 D+00 )
      PARAMETER ( FINCFR = 1.0 D+00 )
      PARAMETER ( ANUC00 = 1.4 D+00 )
      PARAMETER ( RAOB00 = 3.0 D+00 )
      PARAMETER ( ACPAR0 = 0.5D+00 - 0.5D+00 * ( ANUC00 - 1.D+00 )
     &                   / RAOB00 )
*
      COMMON / FKNUCO / HELP (2), HHLP (2), FTVTH (2), FINCX (2),
     &                  EKPOLD (2), BBOLD, ZZOLD, SQROLD, ASEASQ,
     &                  FSPRED, FEX0RD
      DIMENSION KPA (39)
      DIMENSION FRA(2,110)
      REAL RNDM(3)
      SAVE BBONE, FRA, KPA, XIXIXI, ANCOSQ
      LOGICAL LNUCNW, LUFFA
      DATA BBONE / 1.D+00 /
      DATA KPA/1,1,3,3,3,3,3,2,2,3,3,3,3,3,3,3,2,2,3,1,1,2,3,3,3,3,
     &         3,3,3,3,1,2,1,2,2,1,1,1,1/
*  Reduction factors for intran. cascade energy, taken from Alsmiller
*  Incoming baryons
      DATA (FRA(1,I), I=1,107) / .048D0,.076D0,0.10D0,.12D0,
     * .14D0,.16D0,.17D0,.19D0,.21D0,.22D0,.24D0,
     1 .25D0,.26D0,.28D0,.29D0,.30D0,.315D0,.33D0,.34D0,
     * .35D0,.36D0,.38D0,.39D0,.40D0,.41D0,.42D0,
     2 .43D0,.44D0,.46D0,.47D0,.48D0,.49D0,.50D0,.51D0,
     * .52D0,.53D0,.54D0,.55D0,.56D0,.57D0,.58D0,.59D0,
     3 .60D0,.61D0,.62D0,.63D0,.635D0,.64D0,.65D0,.66D0,
     * .67D0,.675D0,.68D0,.69D0,.70D0,.71D0,.715D0,
     4 .72D0,.725D0,.73D0,.735D0,.74D0,.75D0,.755D0,
     * .76D0,.767D0,.77D0,.77D0,.775D0,.78D0,.783D0,
     5 .786D0,.79D0,.795D0,.80D0,.805D0,.81D0,.812D0,
     * .815D0,.82D0,.822D0,.824D0,.825D0,.825D0,.83D0,
     6 .832D0,.834D0,.836D0,.838D0,.84D0,.843D0,.836D0,
     * .849D0,.852D0,.855D0,.856D0,.857D0,.858D0,
     7 .859D0,.860D0,.862D0,.864D0,.866D0,.868D0,.870D0,
     * .872D0,.874D0/
*  Incoming mesons
      DATA (FRA(2,I), I=1,63) / .048D0,.076D0,0.10D0,.12D0,
     * .14D0,.16D0,.17D0,.19D0,.21D0,.22D0,.24D0,
     1 .25D0,.262D0,.274D0,.285D0,.295D0,.305D0,.315D0,
     * .327D0,.340D0,.345D0,.350D0,.360D0,
     2 .370D0,.380D0,.385D0,.390D0,.398D0,.406D0,.415D0,
     * .420D0,.425D0,.430D0,.435D0,.440D0,.446D0,
     3 .452D0,.458D0,.464D0,.470D0,.474D0,.478D0,.482D0,
     * .486D0,.490D0,.495D0,.500D0,.505D0,.510D0,
     4 .515D0,.519D0,.523D0,.526D0,.520D0,.533D0,.537D0,
     * .540D0,.543D0,.547D0,.550D0,.555D0,.559D0,
     5 .5625D0 /
*
*   The calculation of number of collisions inside nucleus is moved
*   from subroutine nucevt; nsea is stored in common corinc
*
      IBPROJ = IBAR (KPROJ)
*  Get the "paprop" index
      IJ     = KPTOIP (KPROJ)
      IJJ    = KPA    (IJ)
*  +-------------------------------------------------------------------*
*  |                                        Incoming baryons
      IF ( IBPROJ .NE. 0 ) THEN
*  |  +----------------------------------------------------------------*
*  |  |
         IF ( IBTAR .LE. 107 ) THEN
            FRAC = FRA ( 1, IBTAR )
*  |  |
*  |  +----------------------------------------------------------------*
*  |  |
         ELSE IF (IBTAR .LE. 206) THEN
            FRAC = 0.739D+00 + 0.00126D+00 * BBTAR
*  |  |
*  |  +----------------------------------------------------------------*
*  |  |
         ELSE
            FRAC = 1.D+00
         END IF
*  |  |
*  |  +----------------------------------------------------------------*
*  |  +----------------------------------------------------------------*
*  |  |   This is a very simple patch to slightly increase the cascade
*  |  |   yield at low energies for antibaryons. The factor 0.5 to
*  |  |   take into account very roughly that not all interactions will
*  |  |   result in complete annihilation
         IF ( IBPROJ .LT. 0 ) THEN
            EKEFF = EKE + 0.5D+00 * AMDISC (KPROJ)
*  |  |
*  |  +----------------------------------------------------------------*
*  |  |
         ELSE
            EKEFF = EKE
         END IF
*  |  |
*  |  +----------------------------------------------------------------*
         STRRED = 0.D+00
*  |  +----------------------------------------------------------------*
*  |  |  Apply a 50 % reduction in the intranuclear cascade probability
*  |  |  for each s/sbar quark of the projectile
         DO 100 IQ = 1,3
            STRRED = STRRED + 0.1666666666666667D+00
     &             * ABS ( IQSCHR ( MQUARK (IQ,IJ) ) )
  100    CONTINUE
*  |  |
*  |  +----------------------------------------------------------------*
*  |   Make the strangeness reduction effective only at low energies
*  |   (where secondaries are few ==> the s composition of the projecti-
*  |    le is relevant)
         STRRED = STRRED * AM (KPROJ) / ( EKE + AM (KPROJ) )
         FRAC   = ( 1.D+00 - STRRED ) * FRAC
*  |
*  +-------------------------------------------------------------------*
*  |                                         Incoming mesons
      ELSE
*  |  +----------------------------------------------------------------*
*  |  |
         IF ( IBTAR .LE. 63 ) THEN
            WEIGH1 = 1.D+00 / ( 1.D+00 + ( 60.D+00 / EKE ) )
            FRAC = FRA ( 2, IBTAR ) * ( ( 1.D+00 + 0.1333D+00 * MAX (
     &             0.D+00, BBTAR - 35.D+00 ) / 28.D+00 ) * WEIGH1 +
     &             1.D+00 - WEIGH1 )
*  |  |
*  |  +----------------------------------------------------------------*
*  |  |
         ELSE IF ( IBTAR .LE. 107 ) THEN
            WEIGH1 = 1.D+00 / ( 1.D+00 + ( 60.D+00 / EKE ) )
            FRAC = ( 0.75D+00 + WEIGH1 * 0.1D+00 ) * FRA ( 1, IBTAR )
*  |  |
*  |  +----------------------------------------------------------------*
*  |  |
         ELSE IF ( IBTAR .LE. 206 ) THEN
            WEIGH1 = 1.D+00 / ( 1.D+00 + ( 60.D+00 / EKE )    )
            FRAC = ( 0.6279D+00 + 0.001077D+00 * BBTAR ) * WEIGH1
     &           + ( 1.D+00 - WEIGH1 ) * ( 0.554D+00 + 0.00095D+00
     &           * BBTAR )
*  |  |
*  |  +----------------------------------------------------------------*
*  |  |
         ELSE
            WEIGH1 = 1.D+00 / ( 1.D+00 + ( 60.D+00 / EKE )    )
            FRAC = 0.75D+00 + 0.1D+00 * WEIGH1
         END IF
*  |  |
*  |  +----------------------------------------------------------------*
         EKEFF = EKE
         STRRED = 0.D+00
*  |  +----------------------------------------------------------------*
*  |  |  Apply a 50 % reduction in the intranuclear cascade probability
*  |  |  for each s/sbar quark of the projectile
         DO 150 IQ = 1,2
            STRRED = STRRED + 0.25D+00
     &             * ABS ( IQSCHR ( MQUARK (IQ,IJ) ) )
  150    CONTINUE
*  |  |
*  |  +----------------------------------------------------------------*
*  |   Make the strangeness reduction effective only at low energies
*  |   (where secondaries are few ==> the s composition of the projecti-
*  |    le is relevant)
         STRRED = STRRED * AM (KPROJ) / ( EKE + AM (KPROJ) )
         FRAC   = ( 1.D+00 - STRRED ) * FRAC
      END IF
*  |
*  +-------------------------------------------------------------------*
*  The following is a reduction factor to bring Hannes's parametri-
*  zations for the average energy carried by cascade nucleons in
*  better agreement with experimental data
      FSPRED = FSPRD0 + ( 1.D+00 - FSPRD0 ) * MAX ( 10.D+00 - EKE,
     &         0.D+00 ) / 10.D+00
      FEX0RD = FSPRD0**ABS(IBPROJ)
*  +-------------------------------------------------------------------*
*  | Check whether it is the same target nucleus of the last call
      IF ( BBTAR .NE. BBOLD .OR. ZZTAR .NE. ZZOLD ) THEN
         LNUCNW = .TRUE.
         SQRAMS = SQRT ( BBTAR )
         ATO1O3 = BBTAR**0.3333333333333333D+00
*        ZTO1O3 = ZZTAR**0.3333333333333333D+00
         BBOLD  = BBTAR
         ZZOLD  = ZZTAR
         SQROLD = SQRAMS
         HKAP  = BBTAR**2 / ( ZZTAR**2 + ( BBTAR - ZZTAR )**2 )
         HHLP (1) = ( HKAP * ZZTAR )**0.3333333333333333D+00 / ATO1O3
         HHLP (2) = ( HKAP * ( BBTAR - ZZTAR ) )
     &              **0.3333333333333333D+00 / ATO1O3
         RDSNUC = R0NUCL * ATO1O3
         RDSCOU = RCCOUL * ATO1O3
         FLKCOU = DOST ( 1, ZZTAR )
*  |  Coulomb barrier "a la Evap"
         VEFFNU (1) = FLKCOU * COULPR * ZZTAR / RDSCOU
         VEFFNU (2) = 0.D+00
         AMRCAV = MAX ( 1.D+00, ( BBTAR - 2.D+00 ) ) * AMUC12
         AMRCSQ = AMRCAV * AMRCAV
*  |  +----------------------------------------------------------------*
*  |  |
         DO 3000 I = 1, 2
            PFRMMX (I) = HHLP (I) * APFRMX
            P2HELP     = PFRMMX (I)**2
            EFRMMX (I) = SQRT ( P2HELP + AMNUSQ (I) ) - AMNUCL (I)
            EFRMAV (I) = 0.3D+00 * P2HELP / AMNUCL (I) * ( 1.D+00
     &                 - P2HELP / ( 5.6D+00 * AMNUSQ (I) ) )
            ERCLAV (I) = 0.3D+00 * P2HELP / AMRCAV * ( 1.D+00
     &                 - P2HELP / ( 5.6D+00 * AMRCSQ ) )
            V0WELL (I) = EFRMMX (I) + EBNDNG (I)
            VEFFNU (I) = VEFFNU (I) + V0WELL (I)
            ERCLMX     = 0.5D+00 * P2HELP / AMRCAV * ( 1.D+00
     &                 - 0.25D+00 * P2HELP / AMRCSQ )
            EKMNNU (I) = VEFFNU (I) + EBNDNG (I)
     &                 - EFRMMX (I) + ERCLMX
            EKMXNU (I) = VEFFNU (I) + EBNDNG (I)
            EKMNAV (I) = VEFFNU (I) + EBNDNG (I)
     &                 - EFRMAV (I) + ERCLAV (I)
            ESWELL (I) = EBNDNG (I) + V0WELL (I)
     &                 - EFRMAV (I) + ERCLAV (I)
            FTVTH  (I) = ( EKMNNU (I) / EFRMMX (I) )**3
            FTVTH  (I) = 0.4D+00 * SQRT  ( FTVTH (I) )
            FINCUP (I) = AKEKA ( I, EKEFF, BBTAR ) * FSPRED
3000     CONTINUE
*  |  |
*  |  +----------------------------------------------------------------*
*  |
*  +-------------------------------------------------------------------*
*  | It is the same target nucleus of the last call to Nucevv/Nucriv
      ELSE
         LNUCNW = .FALSE.
         SQRAMS = SQROLD
         FINCUP (1) = AKEKA ( 1, EKEFF, BBTAR ) * FSPRED
         FINCUP (2) = AKEKA ( 2, EKEFF, BBTAR ) * FSPRED
      END IF
*  |
*  +-------------------------------------------------------------------*
      CALL NIZL (IJ, BBTAR, EKE, PPROJ, SIHA, ZLA)
      CALL NIZL (IJ, BBONE, EKE, PPROJ, SIHN, ZLP)
      ANUAV = BBTAR  * SIHN / SIHA
*  Anuav= average number of collisions in nucleus
      ANUSEA = ANUAV - 1.D+00
*  +-------------------------------------------------------------------*
*  |                Call the function sampling the distribution for the
*  |                number of high energy collisions
      IF ( ANUSEA .GT. 0.D+00 .AND. PPROJ .GT. 5.D+00 ) THEN
         EXPLAM = EXP ( -ANUSEA )
         NSEA   = NUDISV ( ANUAV, IBPROJ, EXPLAM, ASEASQ, APOWER,
     &                     PRZERO ) - 1
         NSEA   = MIN ( NSEA, IBTAR - 1 )
         LUFFA  = .FALSE.
*  |
*  +-------------------------------------------------------------------*
*  |
      ELSE
         ASEASQ = 0.D+00
         ANUSEA = 0.D+00
         APOWER = 1.D+00
         ANUAV  = 1.D+00
         PRZERO = 1.D+00
         ANUCSQ = 1.D+00
         NSEA   = 0
         RATOLD = 0.D+00
         ACFACT = 0.D+00
         BCFACT = 2.D+00
         ANUC11 = ANUC00
         ANUCOR = ANUC00
         ACPARM = 0.D+00
         RRPCO  = 0.D+00
         LUFFA  = .TRUE.
      END IF
*  |
*  +-------------------------------------------------------------------*
      NSEALD = NSEA
*  +-------------------------------------------------------------------*
*  |  Check if the parameterized intranuclear cascade is requested
      IF ( LINCTV ) THEN
         IF ( LUFFA ) GO TO 3500
         CALL GRNDM(RNDM,1)
         RRPCR = RNDM (1)
*  |  +----------------------------------------------------------------*
*  |  |              Check for negative <n_sea^2>
         IF ( ASEASQ .LT. 0.D+00 ) THEN
            ANUCOR = ANUC00 / ( 1.D+00 + 2.D+00 / PPROJ )
            ASEASQ = ANUCOR * ANUAV**2 - 2.D+00 * ANUAV + 1.D+00
            ANUCSQ = ANUCOR * ANUAV**2
            ACPARM = ACPAR0
            ANUC11 = ANUC00
            PRZERO = EXPLAM
*  |  |
*  |  +----------------------------------------------------------------*
*  |  |
         ELSE
            ANUCSQ = ASEASQ + 2.D+00 * ANUAV - 1.D+00
            ANUC11 = ANUCSQ / ANUAV**2
*  |  |  +-------------------------------------------------------------*
*  |  |  |  Check if a power different from 2 must be used for the
*  |  |  |  cascade particle distributions
*  |  |  |  Power = 2:
            IF ( .NOT. LPOWER ) THEN
               APOWER = ANUCSQ
*  |  |  |
*  |  |  +-------------------------------------------------------------*
*  |  |  |
            ELSE
               IPOWER = NINT (-DPOWER)
               IF ( IPOWER .EQ. 17 ) THEN
                  FEX0RD = FEX0RD * FSPRD0
               ELSE IF ( IPOWER .EQ. 8 .OR. IPOWER .EQ. 1 .OR. IPOWER
     &            .GE. 12 ) THEN
                  FEX0RD = FEX0RD * FSPRD0
               END IF
            END IF
*  |  |  |
*  |  |  +-------------------------------------------------------------*
            ANUC11 = MIN ( ANUC11, ANUC00 )
            ANUCOR = ANUC11 / ( 1.D+00 + 5.D+00 * ( ANUC11 - 1.D+00 )
     &             / PPROJ )
            ANUCOR = MAX ( ANUCOR, ONEONE )
            BBCOF  = ( ANUC11 * ANUAV**2 - 2.D+00 * ANUCOR * ANUAV**2
     &             + 2.D+00 * ANUCOR * ANUAV - 1.D+00 )
            BBCOF  = ABS (BBCOF)
            AACOF  = ANUCOR * ( ANUAV - 1.D+00 )**2
            CCCOF  = - ANUAV**2 * ( ANUC11 - ANUCOR )
            IF ( AACOF - CCCOF .GT. ANGLGB ) THEN
               RRPCO = 0.5D+00 * ( - BBCOF + SQRT ( BBCOF**2 - 4.D+00
     &               * AACOF * CCCOF ) ) / AACOF
            ELSE
               RRPCO = 0.D+00
            END IF
            RRPCO  = MIN ( ONEONE, RRPCO )
            RRPCO  = MAX ( ZERZER, RRPCO )
            RRPCO  = 1.D+00 - RRPCO
         END IF
*  |  |
*  |  +----------------------------------------------------------------*
*  |  Supply the fraction of the total kinetic energy to be
*  |  used for intranuclear cascade nucleons
 3500    CONTINUE
         EKUPNU (1) = AINEL (IJJ, 6, EKEFF, BBTAR, SQRAMS) * EKEFF
     &              * FSPRED
         EKUPNU (2) = AINEL (IJJ, 7, EKEFF, BBTAR, SQRAMS) * EKEFF
     &              * FSPRED
         EINCP  = FRAC * EKUPNU (1)
         EINCN  = FRAC * EKUPNU (2)
         FRAMAX = FRAC
*  |  +----------------------------------------------------------------*
*  |  |  Now start Fincup (i) calculation!!!!
         DO 4000 I = 1, 2
            EKPOLD (I) = EKUPNU (I)
            EKUPNU (I) = MAX ( FRAMAX  * EKUPNU (I), EKMXNU (I)/TWOTHI
     &                         , FOUFOU * FINCUP (I) )
            EKUPNU (I) = MIN ( EKUPNU (I), FOUFOU * FINCUP (I) )
            AHELP = - ( EKUPNU (I) - EKMNAV (I) ) / FINCUP (I)
            CHELP = EKMNAV (I) / FINCUP (I)
            DHELP = EKUPNU (I) / FINCUP (I)
            FINC0 = FINCFR + ESWELL (I) / FINCUP (I)
            BHELP = EXP ( AHELP / FINC0 )
            FINCA = FINC0
            FUNCA = - ( CHELP - DHELP * BHELP ) / ( 1.D+00 - BHELP )
            FINCB = 2.D+00 * FINC0
            BHELP = EXP ( AHELP / FINCB )
            FUNCB = FINC0 - FINCB - ( CHELP - DHELP * BHELP ) /
     &            ( 1.D+00 - BHELP )
            ICOU = 1
            FINCLD = FINC0
*  |  |  +-------------------------------------------------------------*
*  |  |  |
3800        CONTINUE
               FINCX (I) = FINCA - FUNCA * ( FINCB - FINCA ) /
     &                   ( FUNCB - FUNCA )
*  |  |  |  +----------------------------------------------------------*
*  |  |  |  |
               IF ( ABS ( FINCX (I) - FINCLD ) .GT. 0.03D+00 .AND.
     &              ICOU .LT. 20 ) THEN
                  ICOU   = ICOU + 1
                  FINCLD = FINCX (I)
*  |  |  |  |  +-------------------------------------------------------*
*  |  |  |  |  |
                  IF ( ABS (FUNCA) .LT. ABS (FUNCB) ) THEN
                     FINCB = FINCLD
                     BHELP = EXP ( AHELP / FINCB )
                     FUNCB = FINC0 - FINCB - ( CHELP - DHELP * BHELP
     &                     ) / ( 1.D+00 - BHELP )
*  |  |  |  |  |
*  |  |  |  |  +-------------------------------------------------------*
*  |  |  |  |  |
                  ELSE
                     FINCA = FINCLD
                     BHELP = EXP ( AHELP / FINCA )
                     FUNCA = FINC0 - FINCA - ( CHELP - DHELP * BHELP
     &                     ) / ( 1.D+00 - BHELP )
                  END IF
*  |  |  |  |  |
*  |  |  |  |  +-------------------------------------------------------*
                  GO TO 3800
*  |  |  |  |
*  |  |  +-<|--<--<--<--<--<
               END IF
*  |  |  |  |
*  |  |  |  +----------------------------------------------------------*
*  |  |  |
*  |  |  +-------------------------------------------------------------*
*  |  |  |
*  |  |  +-------------------------------------------------------------*
            ESLOPE (I) = FINCUP (I) * FINCX (I)
            AHELP      = EXP ( - EKPOLD (I) / FINCUP (I) )
            EKNOLD     = FINCUP (I) - EKPOLD (I) * AHELP /
     &                 ( 1.D+00 - AHELP )
            EXMNAV (I) = EXP ( - EKMNAV (I) / ESLOPE (I) )
            EXMNNU (I) = EXP ( - EKMNNU (I) / ESLOPE (I) )
            EXUPNU (I) = EXP ( - EKUPNU (I) / ESLOPE (I) )
            EKINAV (I) = ESLOPE (I) + ( EKMNAV (I) * EXMNAV (I) -
     &                   EKUPNU (I) * EXUPNU (I) ) / ( EXMNAV (I) -
     &                   EXUPNU (I) )
            EKPOLD (I) = EKPOLD (I) * FRAC
            FINCUP (I) = EKINAV (I) / EKNOLD
4000     CONTINUE
*  |  |
*  |  +----------------------------------------------------------------*
         EINCP = EINCP * FINCUP (1)
         EINCN = EINCN * FINCUP (2)
         AGREYP = EINCP / EKINAV (1)
         AGREYN = EINCN / EKINAV (2)
         AGREYT = AGREYP + AGREYN
         CALL GRNDM(RNDM,1)
         RNPCO  = RNDM ( 1 )
*  |  +----------------------------------------------------------------*
*  |  |  Check if we have to sample from a distribution with <n_casc>
*  |  |  prop. <nu^power> or to <nu>
         IF ( RNPCO .LT. RRPCO ) THEN
*  |  |  +-------------------------------------------------------------*
*  |  |  |  Power .ne. 2:
            IF ( LPOWER ) THEN
*  |  |  |  +----------------------------------------------------------*
*  |  |  |  |  1**Y = 1
               IF ( NSEA .EQ. 0 ) THEN
                  IF ( IPOWER .LT. 7 ) THEN
                     ANCOSQ = 1.D+00
                  ELSE
                     ANCOSQ = FPOWER ( IPOWER, 1, ANUAV )
                  END IF
*  |  |  |  |
*  |  |  |  +----------------------------------------------------------*
*  |  |  |  |  The exponent has a fixed value
               ELSE IF ( DPOWER .GT. 0.D+00 ) THEN
                  ANCOSQ = ( NSEA + 1 )**DPOWER
*  |  |  |  |
*  |  |  |  +----------------------------------------------------------*
*  |  |  |  |  The function Fpower supplies the exponent
               ELSE
                  IF ( IPOWER .LT. 7 ) THEN
                     ANCOSQ = ( NSEA + 1 )**FPOWER ( IPOWER,
     &                          NSEA + 1, ANUAV )
                  ELSE IF ( IPOWER .LT. 11 ) THEN
                     ANCOSQ = ( NSEA + 1 ) * FPOWER ( IPOWER,
     &                        NSEA + 1, ANUAV )
                  ELSE
                     ANCOSQ = ( NSEA + 1 )**FPOWER ( IPOWER,
     &                        NSEA + 1, ANUAV )
                  END IF
               END IF
*  |  |  |  |
*  |  |  |  +----------------------------------------------------------*
*  |  |  |
*  |  |  +-------------------------------------------------------------*
*  |  |  |  Power = 2:
            ELSE
               ANCOSQ = ( NSEA + 1 )**2
            END IF
*  |  |  |
*  |  |  +-------------------------------------------------------------*
            XIXIXI = AGREYP / ( AGREYP + APOWER )
*  |  |
*  |  +----------------------------------------------------------------*
*  |  |
         ELSE
            ANCOSQ = NSEA + 1
            XIXIXI = AGREYP / ( AGREYP + ANUAV )
         END IF
*  |  |
*  |  +----------------------------------------------------------------*
         RRPC0  = ( 1.D+00 - XIXIXI )**ANCOSQ
         RRPCR  = RRPC0
         CALL GRNDM(RNDM,1)
         RNPCR  = RNDM (1)
         IF ( RRPCR .GE. RNPCR ) THEN
            NGREYP = 0
         ELSE
            DO 4600 I = 1, ICHTAR
               RRPC0 = RRPC0 * ( I - 1 + ANCOSQ ) * XIXIXI
     &               / I
               RRPCR = RRPCR + RRPC0
               IF ( RNPCR .LE. RRPCR ) GO TO 4700
 4600       CONTINUE
            I = I - 1
            IF ( BBTAR .GT. 12.D+00 )
     &      WRITE (LUNERR,*)' *** RRPCR,I,NSEA,ANCOSQ*XIXIXI ***',
     &                            RRPCR,I,NSEA,ANCOSQ*XIXIXI
 4700       CONTINUE
            NGREYP = I
         END IF
*  |  +----------------------------------------------------------------*
*  |  |  Check if we have to sample from a distribution with <n_casc>
*  |  |  prop. <nu^power> or to <nu>
         IF ( RNPCO .LT. RRPCO ) THEN
            XIXIXI = AGREYN / ( AGREYN + APOWER )
*  |  |
*  |  +----------------------------------------------------------------*
*  |  |
         ELSE
            ANCOSQ = NSEA + 1
            XIXIXI = AGREYN / ( AGREYN + ANUAV )
         END IF
*  |  |
*  |  +----------------------------------------------------------------*
         RRNC0  = ( 1.D+00 - XIXIXI )**ANCOSQ
         RRNCR  = RRNC0
*  |  Correlate cascade neutrons to cascade protons:
*  |  No correlation
         CALL GRNDM(RNDM,1)
         RNNCR  = RNDM (1)
         RN1GSC = RNPCR
         RN2GSC = RNNCR
         IF ( RRNCR .GE. RNNCR ) THEN
            NGREYN = 0
         ELSE
            DO 4650 I = 1, IBTAR - ICHTAR
               RRNC0 = RRNC0 * ( I - 1 + ANCOSQ ) * XIXIXI
     &               / I
               RRNCR = RRNCR + RRNC0
               IF ( RNNCR .LE. RRNCR ) GO TO 4750
 4650       CONTINUE
            I = I - 1
            IF ( BBTAR .GT. 12.D+00 )
     &      WRITE (LUNERR,*)' *** RRNCR,I,NSEA,ANCOSQ*XIXIXI ***',
     &                            RRNCR,I,NSEA,ANCOSQ*XIXIXI
 4750       CONTINUE
            NGREYN = I
         END IF
         NGREYT = NGREYN + NGREYP
         NGREYN = 0
         NGREYP = 0
         PROBPR = AGREYP / ( AGREYP + AGREYN )
         DO 4760 I = 1, NGREYT
            CALL GRNDM(RNDM,1)
            RNDPPR = RNDM (1)
            IF ( RNDPPR .LT. PROBPR .AND. NGREYP .LT. ICHTAR )
     &         THEN
               NGREYP = NGREYP + 1
            ELSE IF ( NGREYN .LT. IBTAR - ICHTAR ) THEN
               NGREYN = NGREYN + 1
            ELSE
               NGREYP = NGREYP + 1
            END IF
 4760    CONTINUE
         FRAMAX = FRAC
*  |  The number of grey particles is now correlated to the number
*  |  of high energy collisions
         EINCP = NGREYP * EKINAV (1)
         EINCN = NGREYN * EKINAV (2)
         TVTENT = ( NSEA + 1 ) * ( AV0WEL - AEFRMA )
*  |
*  +-------------------------------------------------------------------*
*  |
      ELSE
         EINCP = 0.D+00
         EINCN = 0.D+00
         EKUPNU (1) = 1.D+01
         EKUPNU (2) = 1.D+01
         ESLOPE (1) = 0.D+00
         ESLOPE (2) = 0.D+00
         EKINAV (1) = 1.D+10
         EKINAV (2) = 1.D+10
         FRAC  = 0.D+00
         PCR   = 1.D+00
         FRAINC = 0.D+00
         TVTENT = 0.D+00
      END IF
*  |
*  +-------------------------------------------------------------------*
      PCROLD = PCR
      EKTRIA = EKE - EINCP - EINCN - TVTENT
      IF ( EKTRIA .LE. 0.5D+00 * EKE ) THEN
         EKTRIA = 0.5D+00 * EKE
      END IF
      ETRIAL = EKTRIA + AM (KPROJ)
*  +-------------------------------------------------------------------*
*  |
      IF ( ETRIAL .LT. ECUTRF ) THEN
         PTRIAL = SQRT ( ETRIAL**2 - AM (KPROJ)**2 )
         XCUTFF =  0.3D+00 * ETHSEA / EKTRIA
*        XCUTFF =  0.3D+00 / PTRIAL
*  |
*  +-------------------------------------------------------------------*
*  |
      ELSE
         PTRIAL = ETRIAL
         UMO = SQRT ( 2.D+00*AM(1) * ETRIAL + AM (KPROJ)**2 + AM(1)**2 )
         XCUTFF = 0.30D+00 * ETHSEA * UMO / ( UMOREF * EKTRIA )
      END IF
*  |
*  +-------------------------------------------------------------------*
 200  CONTINUE
      IF ( NSEA .EQ. 0 ) GO TO 1000
* *** Sample X-values of nsea sea-quark-antiquark pairs
      NC3 = 0
*  +-------------------------------------------------------------------*
*  |
 300  CONTINUE
* *** Sea distribution X**(-1)*(1-X)**5
         NC3 = NC3 + 1
*  |  +----------------------------------------------------------------*
*  |  |
         IF ( NC3 .GT. 10 ) THEN
            NSEA = NSEA - 1
            GO TO 200
         END IF
*  |  |
*  |  +----------------------------------------------------------------*
         XO = 1.D+00
         UNOSEA = 2.0D+00
         GAMSEA = 0.05D+00
         XSEAMX = 1.0D+00
*  |  +----------------------------------------------------------------*
*  |  |
         DO 400 I=1,NSEA
            NCOU= 0
  22        CONTINUE
            NCOU = NCOU + 1
*  |  |  +-------------------------------------------------------------*
*  |  |  |
            IF ( NCOU .LT. 11 ) THEN
               XSEA(I) = BETRST ( GAMSEA, UNOSEA, XCUTFF, XSEAMX )
               IF ( XSEA(I) .LT. XCUTFF ) GO TO 22
            END IF
*  |  |  |
*  |  |  +-------------------------------------------------------------*
            NCOU = 0
  23        CONTINUE
            NCOU = NCOU + 1
*  |  |  +-------------------------------------------------------------*
*  |  |  |
            IF ( NCOU .LT. 11 ) THEN
               XASEA(I) = BETRST ( GAMSEA, UNOSEA, XCUTFF, XSEAMX )
               IF ( XASEA(I) .LT. XCUTFF ) GO TO 23
            END IF
*  |  |  |
*  |  |  +-------------------------------------------------------------*
            XO = XO * ( 1.D+00 - XSEA (I) - XASEA(I) )
 400     CONTINUE
*  |  |
*  |  +----------------------------------------------------------------*
*  |  +----------------------------------------------------------------*
*  |  |
         IF ( XO * PTRIAL .LT. 4.D+00 ) THEN
            NSEA = NSEA - 1
            IF ( NSEA  .LE. 0 ) GO TO 1000
            GO TO 300
*  |  |
*  +--+--<--<--<--<--< go to resample
*  |  |
         END IF
*  |  |
*  |  +----------------------------------------------------------------*
         ILOW = 0
*  |  +----------------------------------------------------------------*
*  |  |
         DO 500 I = 1, NSEA
            EMSEA = EKTRIA * ( XSEA (I) + XASEA (I) )
*  |  |  +-------------------------------------------------------------*
*  |  |  |
            IF ( EMSEA .LE. ETHSEA + AM (23) ) THEN
               ILOW = ILOW + 1
*  |  |  |  +----------------------------------------------------------*
*  |  |  |  |
               IF ( I .LT. NSEA ) THEN
                  II = I + 1
                  XSEA  (II) = XSEA  (II) + XSEA  (I)
                  XASEA (II) = XASEA (II) + XASEA (I)
*  |  |  |  |
*  |  |  |  +----------------------------------------------------------*
*  |  |  |  |
               ELSE IF ( I - ILOW .GT. 0 ) THEN
                  II = I - ILOW
                  EKTRIA = EKTRIA * ( 1.D+00 - XSEA(I) - XASEA (I) )
                  XSEA  (II) = XSEA  (II) + XSEA  (I)
                  XASEA (II) = XASEA (II) + XASEA (I)
               END IF
*  |  |  |  |
*  |  |  |  +----------------------------------------------------------*
*  |  |  |
*  |  |  +-------------------------------------------------------------*
*  |  |  |
            ELSE
               XSEA  (I-ILOW) = XSEA  (I)
               XASEA (I-ILOW) = XASEA (I)
               EKTRIA = EKTRIA * ( 1.D+00 - XSEA(I) - XASEA (I) )
            END IF
*  |  |  |
*  |  |  +-------------------------------------------------------------*
 500     CONTINUE
*  |  |
*  |  +----------------------------------------------------------------*
         NSEA = NSEA - ILOW
*  |
*  +-------------------------------------------------------------------*
1000  CONTINUE
*  +-------------------------------------------------------------------*
*  |     Now sample the target nucleons for the high energy collisions
      DO 1200 I = 1, NSEA + 1
         ZAPU = ZNOW / ANOW
*  |  +----------------------------------------------------------------*
*  |  |  Wounded nucleon selection:
*  |  |  Kt is the index of the target nucleon (1=proton, 8=neutron)
         CALL GRNDM(RNDM,1)
         IF ( RNDM(1) .LE. ZAPU ) THEN
            IJTARG (I) = 1
            ZNOW  = ZNOW - 1.D0
            KTARP = KTARP + 1
*  |  |
*  |  +----------------------------------------------------------------*
*  |  |
         ELSE
            IJTARG (I) = 8
            KTARN = KTARN + 1
         END IF
*  |  |
*  |  +----------------------------------------------------------------*
         ANOW   = ANOW - 1.D0
1200  CONTINUE
*  |
*  +-------------------------------------------------------------------*
      ZNCOLL = KTARP
      ANCOLL = NSEA + 1
      AHELP  = EKMNNU (1) / ESLOPE (1)
      AHELP  = AHELP * FTVTH (1) * ( 1.D+00 + 0.3333333333333333D+00
     &       * AHELP ) / ( 1.D+00 - EXUPNU (1) )
      AHELPP = AHELP
      AHELP  = EKMNNU (2) / ESLOPE (2)
      AHELP  = AHELP * FTVTH (2) * ( 1.D+00 + 0.3333333333333333D+00
     &       * AHELP ) / ( 1.D+00 - EXUPNU (2) )
      AHELPN = AHELP
      FEXTRA = 0.3D+00 * AGREYP
      TVCHCP = EFRMMX (1) - EFRMAV (1) + EBNDNG (1)
      TVGRYP = AHELPP * ( EINCP + FEXTRA * EKINAV (1) ) * AGREYP
     &       / ( AGREYP + FEXTRA ) * FEX0RD
      NHELP  = INT ( TVGRYP / TVCHCP )
      TVGRE0 = NHELP * TVCHCP
      PROB0  = ( TVGRYP - TVGRE0 ) / TVCHCP
      CALL GRNDM(RNDM,1)
      IF ( RNDM (1) .LT. PROB0 ) THEN
         CALL GRNDM(RNDM,3)
         P2HELP = ( PFRMMX (1) * MAX ( RNDM (1),
     &              RNDM (2), RNDM (3) ) )**2
         TVGRYP = 0.5D+00 * P2HELP / AMNUCL (1) * ( 1.D+00
     &          - 0.25D+00 * P2HELP / AMNUSQ (1) )
         TVGRYP = EFRMMX (1) - TVGRYP + EBNDNG (1) + TVGRE0
      ELSE
         TVGRYP = TVGRE0
      END IF
      FEXTRA = 0.3D+00 * AGREYN
      TVCHCN = EFRMMX (2) - EFRMAV (2) + EBNDNG (2)
      TVGRYN = AHELPN * ( EINCN + FEXTRA * EKINAV (2) ) * AGREYN
     &       / ( AGREYN + FEXTRA ) * FEX0RD
      NHELP  = INT ( TVGRYN / TVCHCN )
      TVGRE0 = NHELP * TVCHCN
      PROB0  = ( TVGRYN - TVGRE0 ) / TVCHCN
      CALL GRNDM(RNDM,1)
      IF ( RNDM (1) .LT. PROB0 ) THEN
         CALL GRNDM(RNDM,3)
         P2HELP = ( PFRMMX (2) * MAX ( RNDM (1),
     &              RNDM (2), RNDM (3) ) )**2
         TVGRYN = 0.5D+00 * P2HELP / AMNUCL (2) * ( 1.D+00
     &          - 0.25D+00 * P2HELP / AMNUSQ (2) )
         TVGRYN = EFRMMX (2) - TVGRYN + EBNDNG (2) + TVGRE0
      ELSE
         TVGRYN = TVGRE0
      END IF
      TVGRE0 = TVGRYP + TVGRYN
      TVGREY = 0.D+00
      NDIFFT = NGREYT - NINT (ANOW)
      IF ( NINT (ZNOW) - NGREYP .LT. 0 ) THEN
         NDIFFP = NGREYP - NINT ( ZNOW )
         NGREYP = NGREYP - NDIFFP
         EINCP  = EINCP - NDIFFP * EKINAV (1)
         NGREYN = NGREYN + NDIFFP
         EINCN  = EINCN + NDIFFP * EKINAV (2)
         IF ( NINT (ANOW-ZNOW) - NGREYN .LT. 0 ) THEN
            NDIFFN = NGREYN - NINT ( ANOW - ZNOW )
            NGREYN = NGREYN - NDIFFN
            EINCN  = EINCN - NDIFFN * EKINAV (2)
         END IF
      ELSE IF ( NINT (ANOW-ZNOW) - NGREYN .LT. 0 ) THEN
         NDIFFN = NGREYN - NINT ( ANOW - ZNOW )
         NGREYN = NGREYN - NDIFFN
         EINCN  = EINCN - NDIFFN * EKINAV (2)
         NGREYP = NGREYP + NDIFFN
         EINCP  = EINCP + NDIFFN * EKINAV (1)
         IF ( NINT (ZNOW) - NGREYP .LT. 0 ) THEN
            NDIFFP = NGREYP - NINT ( ZNOW )
            NGREYP = NGREYP - NDIFFP
            EINCP  = EINCP - NDIFFP * EKINAV (1)
         END IF
      END IF
      NGREYT = NGREYP + NGREYN
      RETURN
      END
