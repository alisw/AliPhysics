*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:19:54  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.42  by  S.Giani
*-- Author :
*$ CREATE CORRIN.FOR
*COPY CORRIN
*                                                                      *
*=== corrin ===========================================================*
*                                                                      *
      SUBROUTINE CORRIN ( ZZTAR, BBTAR, KPROJ, PPROJ, EKE )
 
#include "geant321/dblprc.inc"
#include "geant321/dimpar.inc"
#include "geant321/iounit.inc"
*
*----------------------------------------------------------------------*
*                                                                      *
*     Created on  10  may  1990    by    Alfredo Ferrari & Paola Sala  *
*                                                   Infn - Milan       *
*                                                                      *
*     Last change on 23-mar-93     by    Alfredo Ferrari               *
*                                                                      *
*    This version has been developed by A. Ferrari starting from the   *
*    the one by J. M. Zazula: it includes a few differences mainly in  *
*    the correlation since now we try to get the excitation energy     *
*    from a correct energy balance                                     *
*----------------------------------------------------------------------*
*
#include "geant321/balanc.inc"
#include "geant321/corinc.inc"
#include "geant321/nucdat.inc"
#include "geant321/parevt.inc"
#include "geant321/part.inc"
#include "geant321/resnuc.inc"
      PARAMETER ( ALFMAX = 5.00D+00 )
      PARAMETER ( FINCFR = 1.00D+00 )
      PARAMETER ( TMNOLD = 0.02D+00 )
      COMMON / FKNUCO / HELP (2), HHLP (2), FTVTH (2), FINCX (2),
     &                  EKPOLD (2), BBOLD, ZZOLD, SQROLD, ASEASQ,
     &                  FSPRED, FEX0RD
      DIMENSION KPA (26)
      DIMENSION FRA(2,110), ESLOLD (2), AVMULT (2), EXDUMM (2),
     &          V0WSAV (2), VEWSAV (2)
      REAL RNDM(3)
      LOGICAL LSTOP, LNUCNW, LLLIN, LLPOW
      SAVE FRA, FRAC, FRAC0, FRAMAX, KPA, IJJ, LSTOP, V0WSAV, VEWSAV
      DATA KPA/1,1,3,3,3,3,3,2,2,3,3,3,3,3,3,3,2,2,3,1,1,2,3,3,3,3/
*
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
      IBPROJ = IBAR (KPROJ)
      IJJ    = KPA  (KPROJ)
      ANUAV  = 1.D+00
      ANCOLL = 1.D+00
      ANUSEA = 1.D+00
      NSEA   = 0
      LSTOP = .FALSE.
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
*  |
*  +-------------------------------------------------------------------*
*  |                                         Incoming mesons
      ELSE
*  |  +----------------------------------------------------------------*
*  |  |
         IF ( IBTAR .LE. 63 ) THEN
            FRAC = FRA ( 2, IBTAR ) * ( 1.D+00 + 0.1333D+00 * MAX (
     &             0.D+00, BBTAR - 35.D+00 ) / 28.D+00 )
*  |  |
*  |  +----------------------------------------------------------------*
*  |  |
         ELSE IF ( IBTAR .LE. 107 ) THEN
            FRAC = 0.85D+00 * FRA ( 1, IBTAR )
*  |  |
*  |  +----------------------------------------------------------------*
*  |  |
         ELSE IF ( IBTAR .LE. 206 ) THEN
            FRAC = 0.6279D+00 + 0.001077D+00 * BBTAR
*  |  |
*  |  +----------------------------------------------------------------*
*  |  |
         ELSE
            FRAC = 0.85D+00
         END IF
*  |  |
*  |  +----------------------------------------------------------------*
      END IF
*  |
*  +-------------------------------------------------------------------*
*  +-------------------------------------------------------------------*
*  |
      IF ( BBTAR .NE. BBOLD .OR. ZZTAR .NE. ZZOLD ) THEN
         LNUCNW = .TRUE.
*  |  Supply the fraction of the total kinetic energy to be
*  |  used for intranuclear cascade nucleons
         SQRAMS = SQRT ( BBTAR )
         ATO1O3 = BBTAR**0.3333333333333333D+00
*        ZTO1O3 = ZZTAR**0.3333333333333333D+00
         HKAP   = BBTAR**2 / ( ZZTAR**2 + ( BBTAR - ZZTAR )**2 )
         HHLP (1) = ( HKAP * ZZTAR )**0.3333333333333333D+00 / ATO1O3
         HHLP (2) = ( HKAP * ( BBTAR - ZZTAR ) )
     &              **0.3333333333333333D+00 / ATO1O3
         RDSNUC = R0NUCL * ATO1O3
         RDSCOU = RCCOUL * ATO1O3
         FLKCOU = DOST ( 1, ZZTAR )
         VEFFNU (1) = COULPR * FLKCOU * ZZTAR / RDSCOU
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
            V0WSAV (I) = V0WELL (I)
            VEWSAV (I) = VEFFNU (I)
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
3000     CONTINUE
*  |  |
*  |  +----------------------------------------------------------------*
         BBOLD  = BBTAR
         ZZOLD  = ZZTAR
         SQROLD = SQRAMS
*  |
*  +-------------------------------------------------------------------*
*  |
      ELSE
         V0WSAV (1) = V0WELL (1)
         VEWSAV (1) = VEFFNU (1)
         V0WSAV (2) = V0WELL (2)
         VEWSAV (2) = VEFFNU (2)
         LNUCNW = .FALSE.
      END IF
*  |
*  +-------------------------------------------------------------------*
      FRAC0 = FRAC
      RETURN
      ENTRY CORSTP ( EKEFF )
         LSTOP = .TRUE.
      ENTRY CORRNC ( EKEFF )
      WEIGH1 = MIN ( 1.D+00, EKEFF / 4.D+00 )
      FRAC   = ( WEIGH1 * 2.5D+00 + ( 1.D+00 - WEIGH1 ) ) * FRAC0
      FRAMAX = 2.5D+00
*  The following is a reduction factor to bring Hannes's parametri-
*  zations for the average energy carried by cascade nucleons in
*  better agreement with experimental data
      FSPRED = FSPRD0 + ( 1.D+00 - FSPRD0 ) * MAX ( 8.D+00 - EKEFF,
     &         0.D+00 ) / 8.D+00
      FSPRED = FSPRED * MIN ( EKEFF / 0.8D+00, 1.D+00 )
      TMPFSP = 1.25D+00 * FSPRD0
      FSPRED = MIN ( FSPRED, TMPFSP )
      AVMULT (1) = FRAC * BNKEKA ( 1, EKEFF, BBOLD, SQROLD )
      AVMULT (2) = FRAC * BNKEKA ( 2, EKEFF, BBOLD, SQROLD )
      EKUPNU (1) = BEKEKA ( 2, EKEFF, BBOLD, SQROLD )
      EKUPNU (2) = BEKEKA ( 3, EKEFF, BBOLD, SQROLD )
      EINCP  = FRAC * EKUPNU (1)
      EINCN  = FRAC * EKUPNU (2)
      EKMAX  = EKEFF - EBNDAV
      EKUPN0 = MAX ( EKUPNU (1), EKUPNU (2) )
*  +-------------------------------------------------------------------*
*  |
      DO 4000 I = 1, 2
         FINCUP (I) = BKEKA ( I, EKEFF, BBOLD )
         RATRAT  =  MIN ( 1.D+00, 0.5D+00 * EKMAX / FINCUP (I) )
         TMPFIN  = 0.5D+00 * EKMAX
         FINCUP (I) = MIN ( FINCUP (I), TMPFIN )
         EKPOLD (I) = EKUPNU (I) * FRAC * RATRAT
         EKUPNU (I) = FINCUP (I)
         TMPEKU     = 1.5D+00 * EKMXNU (I)
         EKUPNU (I) = MAX ( FRAMAX * EKUPNU (I), TMPEKU )
         TMPEKU     = 0.7D+00 * EKEFF
         EKUPNU (I) = MIN ( EKUPNU (I), TMPEKU )
         ESLOLD (I) = MAX ( FINCUP (I), TMNOLD )
         FINCX  (I) = ( ESLOLD (I) + ESWELL (I) ) / ESLOLD (I) * FSPRED
         ESLOPE (I) = MIN ( ESLOLD (I) * FINCX (I), EKMAX )
         AHELP      = EXP ( - EKPOLD (I) / ESLOLD (I) )
         EKNOLD     = ESLOLD (I) - EKPOLD (I) * AHELP /
     &              ( 1.D+00 - AHELP )
         EXMNAV (I) = EXP ( - EKMNAV (I) / ESLOPE (I) )
         EXMNNU (I) = EXP ( - EKMNNU (I) / ESLOPE (I) )
         EXUPNU (I) = EXP ( - EKUPNU (I) / ESLOPE (I) )
         EKINAV (I) = ESLOPE (I) + ( EKMNAV (I) * EXMNAV (I) -
     &                EKUPNU (I) * EXUPNU (I) ) / ( EXMNAV (I) -
     &                EXUPNU (I) )
         FINCUP (I) = AVMULT (I) * EKINAV (I) / EKPOLD (I)
*        FINCUP (I) = EKINAV (I) / EKNOLD
4000  CONTINUE
*  |
*  +-------------------------------------------------------------------*
      EINCP = EKPOLD (1)
      EINCN = EKPOLD (2)
      IF ( LSTOP ) RETURN
*  Power sampling
      LLPOW  = .TRUE.
      LLLIN  = .FALSE.
*  Take into account that a fraction 1/A^2/3 is lost because of
*  peripheral collisions
      PERCOR = ATO1O3 * ATO1O3
      PERCOR = PERCOR / ( PERCOR - 1.D+00 )
      EINCP0 = EINCP  * FINCUP (1) * PERCOR
      EINCN0 = EINCN  * FINCUP (2) * PERCOR
      EINCT  = EINCP0 + EINCN0
* For the other distr,:
      EINCMX = EKEFF - EBNDAV - ( AV0WEL - AEFRMA ) / ( 1.D+00 +
     &         4.D+00 * ( AV0WEL - AEFRMA ) / EKEFF )
      EIUSE  = 0.5D+00 * ( EINCMX + EKMAX )
      AHNORM = 0.D+00
      EINCM0 = - AINFNT
*  +-------------------------------------------------------------------*
*  |  Compute alfa:
      IF ( EINCMX .GT. 2.1D+00 * EINCT ) THEN
         EINFIN = MIN ( 0.95D+00 * EINCMX, ( ALFMAX + 2.D+00 ) * EINCT )
         ALFA   = EINFIN / EINCT - 2.D+00
*  |
*  +-------------------------------------------------------------------*
*  |  Energy interval too small, sample E uniformly
      ELSE
         LLPOW  = .FALSE.
      END IF
*  |
*  +-------------------------------------------------------------------*
      EINCMN = 0.D+00
*  +-------------------------------------------------------------------*
*  | Energy is so small that we cannot produce cascade nucleons,
*  | sample an excitation energy and return
      IF ( EINCMX .LE. EINCMN + 0.25D+00 * EBNDAV ) THEN
         ISAMPL = 50
*  |
*  +-------------------------------------------------------------------*
*  |
      ELSE
         ISAMPL = 0
      END IF
*  |
*  +-------------------------------------------------------------------*
*  +-------------------------------------------------------------------*
*  |  Now make the linear sampling for Eincp, Eincn, on [0,Eicnmx]
*  |  according to P(E)dE=A(Efin-E)^a dE
5000  CONTINUE
         ISAMPL = ISAMPL + 1
         CALL GRNDM(RNDM,1)
         IF ( LLPOW ) THEN
            EIHELP = EINFIN * ( 1.D+00 - RNDM (1)**(1.D+00
     &             /(ALFA+1.D+00)) )
         ELSE IF ( LLLIN ) THEN
            EINCMX = EINCM0
            EIHELP = EINFIN - SQRT ( EINFIN**2 - 2.D+00 * RNDM (1)
     &             / AHNORM )
         ELSE
            EIHELP = EINCMX * RNDM (1)
         END IF
*  |  *** end linear sampling ***
         EINCP  = EINCP0 * EIHELP / EINCT
         EINCN  = EINCN0 * EIHELP / EINCT
         AHELP  = EKMNNU (1) / ESLOPE (1)
         AHELP  = AHELP * FTVTH (1) * EINCP * ( 1.D+00 + AHELP *
     &            0.3333333333333333D+00 ) / ( 1.D+00 - EXUPNU (1) )
         TVGRE0 = AHELP
         AHELP  = EKMNNU (2) / ESLOPE (2)
         AHELP  = AHELP * FTVTH (2) * EINCN * ( 1.D+00 + AHELP *
     &            0.3333333333333333D+00 ) / ( 1.D+00 - EXUPNU (2) )
         TVGRE0 = TVGRE0 + AHELP
*  |  +----------------------------------------------------------------*
*  |  |  Energy is so small that we cannot produce cascade nucleons,
*  |  |  sample an excitation energy and return
         IF ( ISAMPL .GE. 50 ) THEN
            CALL GRNDM(RNDM,3)
            PRNDM  = MAX ( RNDM (1), RNDM (2), RNDM (3) )
            TVGRE0 = MAX ( AV0WEL - PRNDM**2 * AEFRMX - EBNDAV, ZERZER )
            EINCP  = 0.D+00
            EINCN  = 0.D+00
            TVGREY = 0.D+00
            RETURN
         END IF
*  |  |
*  |  +----------------------------------------------------------------*
*  | Changed:
         IF ( EIHELP + TVGRE0 .GE. EIUSE )  GO TO 5000
*  |
*  +--<--<--<--<--<  Resampling!
*  |
*  +-------------------------------------------------------------------*
      PCR    = EIHELP / EINCT
      FRAINC = FRAC * PCR
      TVGREY = 0.D+00
      AGREYP = EINCP / EKINAV (1)
      AGREYN = EINCN / EKINAV (2)
* ==== Discretize the distribution !!! ==== *
      CALL GRNDM(RNDM,2)
      NGREYP = INT ( AGREYP )
      IF ( RNDM(1) .LT. AGREYP - NGREYP ) NGREYP = NGREYP + 1
      NGREYN = INT ( AGREYN )
      IF ( RNDM(2) .LT. AGREYN - NGREYN ) NGREYN = NGREYN + 1
* ====
      AGREYT = AGREYP + AGREYN
      NGREYT = NGREYP + NGREYN
      PPROCS = AGREYP / AGREYT
      NGREYP = 0
      NGREYN = 0
      EINCP  = 0.D+00
      EINCN  = 0.D+00
      IF ( NGREYT .LE. 0 ) GO TO 9000
      ARDUMM = MAX ( 1.D+00, BBOLD - 2.D+00 )
      CALL RBKINI ( 1, .TRUE., EXDUMM, TKDUMM, TSDUMM,
     &               PSDUMM, ARDUMM, TRDUMM )
      V0SAV1 = V0WELL (1)
      V0SAV2 = V0WELL (2)
      VESAV1 = VEFFNU (1)
      VESAV2 = VEFFNU (2)
      V0WELL (1) = V0WSAV (1)
      V0WELL (2) = V0WSAV (2)
      VEFFNU (1) = VEWSAV (1)
      VEFFNU (2) = VEWSAV (2)
      IRETRY = 0
      DO 8000 I = 1, NGREYT
         CALL GRNDM(RNDM,1)
         RNDMPR = RNDM (1)
         ARDUMM = MAX ( 1.D+00, ARDUMM - 1.D+00 )
 7500    CONTINUE
         IF ( RNDMPR .LT. PPROCS ) THEN
            NGREYP = NGREYP + 1
            CALL RBKINI ( 1, .FALSE., EXDUMM, TKDUMM, TSDUMM,
     &                    PSDUMM, ARDUMM, TRDUMM )
            EINCP = EINCP + TKDUMM
            IF ( EINCP + EINCN .GT. EINCMX .AND. IRETRY .LT. 5 ) THEN
               EINCP = EINCP - TKDUMM
               CALL RBKMIN (1)
               IRETRY = IRETRY + 1
               GO TO 7500
            END IF
         ELSE
            NGREYN = NGREYN + 1
            CALL RBKINI ( 2, .FALSE., EXDUMM, TKDUMM, TSDUMM,
     &                    PSDUMM, ARDUMM, TRDUMM )
            EINCN = EINCN + TKDUMM
            IF ( EINCP + EINCN .GT. EINCMX .AND. IRETRY .LT. 5 ) THEN
               EINCN = EINCN - TKDUMM
               CALL RBKMIN (2)
               IRETRY = IRETRY + 1
               GO TO 7500
            END IF
         END IF
 8000 CONTINUE
      V0WELL (1) = V0SAV1
      V0WELL (2) = V0SAV2
      VEFFNU (1) = VESAV1
      VEFFNU (2) = VESAV2
 9000 CONTINUE
      EINCT  = EINCP  + EINCN
      IF ( EINCT + TVGRE0 .GT. EIUSE ) THEN
         TVGRE0 = 0.D+00
         EIUSE = MIN ( ONEONE, EINCMX / EINCT )
         EINCP = EINCP * EIUSE
         EINCN = EINCN * EIUSE
         EINCT = EINCP + EINCN
      END IF
      RETURN
*=== End of subroutine corrin =========================================*
      END
