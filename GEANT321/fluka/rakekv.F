*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:19:58  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.44  by  S.Giani
*-- Author :
*$ CREATE RAKEKV.FOR
*COPY RAKEKV
*                                                                      *
*=== rakekv ===========================================================*
*                                                                      *
      SUBROUTINE RAKEKV ( IT, EXSOP, BBTAR, TKIN, TSTRCK, PSTRCK,
     &                    ARECL, TKRECL, EFERMI, CDE, SDE )
 
#include "geant321/dblprc.inc"
#include "geant321/dimpar.inc"
#include "geant321/iounit.inc"
*
*----------------------------------------------------------------------*
*     version by                     Alfredo Ferrari
*                                    INFN - Milan
*     last change 03 january 93 by   Alfredo Ferrari
*                                    INFN - Milan
*
*     To be called from the high energy production
*
*     this is a subroutine of fluka to sample intranuclear cascade
*     particles: it is based on the old Rakeka from J. Ranft
*
*     input variables:
*        it    = type of the secondary requested; 1=proton, 2=neutron
*        bbtar = atomic weight of the medium
*
*     output variables:
*        tkin  = kinetic energy of the secondary in GeV before applying
*                the nuclear well (and eventually the Coulomb barreer)
*        tstrck= kinetic energy of the secondary in GeV
*        pstrck= momentum of the secondary in GeV/c
*        cde,sde = cosine and sine of the angle between the
*                  directions of the primary
*                  and secondary particles
*
*********************************************************************
*
#include "geant321/nucdat.inc"
#include "geant321/parevt.inc"
      COMMON / FKNUCO / HELP (2), HHLP (2), FTVTH (2), FINCX (2),
     &                  EKPOLD (2), BBOLD, ZZOLD, SQROLD, ASEASQ,
     &                  FSPRED, FEX0RD
      COMMON / FKCASF / PKFRMI, COSTH, PKIN
      DIMENSION EXSOP (2)
      REAL RNDM(4), RNGSS, DUMNOR
*  In this version the low energy component has been suppressed, since
*  they are now produced by the evaporation model, A. Ferrari.
*  The parameters needed for the sampling have been already set in
*  Corrin
      EXSOP (IT) = 0.D+00
*  +-------------------------------------------------------------------*
*  |
 100  CONTINUE
*  |  Sample the Fermi momentum
         CALL GRNDM(RNDM,4)
         PKFRMI = PFRMMX (IT) * MAX ( RNDM (1), RNDM (2),
     &            RNDM (3) )
         PKFRSQ = PKFRMI * PKFRMI
         TKIN   = - ESLOPE (IT) * LOG ( EXMNNU (IT) - RNDM (4) * (
     &              EXMNNU (IT) - EXUPNU (IT) ) )
         TKRECL = 0.5D+00 * PKFRSQ / ( AMUC12 * ARECL ) * ( 1.D+00 -
     &            0.25D+00 * PKFRSQ / ( ARECL**2 * AMUCSQ ) )
         EFERMI = SQRT ( AMNUSQ (IT) + PKFRSQ )
         TSTRCK = EFERMI + TKIN - AMNUCL (IT) - V0WELL (IT) - TKRECL -
     &            EBNDNG (IT)
*  |  +----------------------------------------------------------------*
*  |  |  Record the energy spent without emitting nucleons
         IF ( TSTRCK .LE. VEFFNU (IT) - V0WELL (IT) ) THEN
*  |  |  Reduce that energy according to Fspred:
            EXSOP (IT) = EXSOP (IT) + TKIN * FEX0RD
            GO TO 100
*  +-<|--<--<--<--<--< go to resampling
         END IF
*  |  |
*  |  +----------------------------------------------------------------*
*  |
*  +-------------------------------------------------------------------*
      PSTRCK = SQRT ( TSTRCK * ( TSTRCK + 2.D+00 * AMNUCL (IT) ) )
*   Sample the angle between the incident particle and the Fermi
*   momentum
*   Solution assuming that the momentum transfer is equal to the
*   the momentum loss of the projectile for the energy loss Tkin,
*   roughly given by Tkin again. Use this with the "usual" eventv
      PKIN   = TKIN
      PKIN2  = PKIN * PKIN
      CALL GRNDM(RNDM,1)
      COSTH  = 2.D+00 * RNDM (1) - 1.D+00
      COSDE2 = PKFRMI * COSTH + PKIN
      CSIGN  = SIGN ( ONEONE, COSDE2 )
      COSDE2 = COSDE2 * COSDE2 / ( PKIN2 + PKFRSQ +
     &           2.D+00 * COSTH * PKIN * PKFRMI )
      CDE = CSIGN * SQRT ( COSDE2 )
* Original
      SDE0 = MAX ( ONEONE - CDE, ANGLGB )
      SDE1 = SDE0
      FCORR0 = 0.10D+00
      TMPSDE = 3.5D+00 * SDE0
      FCORR  = ATO1O3 * MIN ( FCORR0, TMPSDE )
 2000 CONTINUE
         CALL GRANOR ( RNGSS, DUMNOR )
         SDE = SDE0 + FCORR * RNGSS
         CDE = 1.D+00 - SDE
      IF ( ABS ( CDE ) .GE. 1.D+00 ) GO TO 2000
      SDE = SQRT ( ( 1.D+00 - CDE ) * ( 1.D+00 + CDE ) )
      RETURN
      END
