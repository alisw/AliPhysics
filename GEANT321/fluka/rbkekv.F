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
*$ CREATE RBKEKV.FOR
*COPY RBKEKV
*                                                                      *
*=== rbkekv ===========================================================*
*                                                                      *
      SUBROUTINE RBKEKV ( IT, EXSOP, TO, AMSS, TKIN, TSTRCK,
     &                    PSTRCK, ARECL, TKRECL, COD, SID )
 
#include "geant321/dblprc.inc"
#include "geant321/dimpar.inc"
#include "geant321/iounit.inc"
*
*----------------------------------------------------------------------*
*     version by                     Alfredo Ferrari
*                                    INFN - Milan
*     last change 19 april 93 by     Alfredo Ferrari
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
*        sid,cod = sine and cosine of the angle between
*                  the directions of the primary
*                  and secondary particles
*
*********************************************************************
*
      PARAMETER ( PI   = PIPIPI )
      PARAMETER ( PIO2 = PI / 2.D+00 )
#include "geant321/nucdat.inc"
      LOGICAL LINIT, LZEROI
      DIMENSION EXSOP (2), AMMOLD (2)
      DIMENSION PFINIT (50,2), TKINIT (50,2), TSINIT (50,2),
     &          TRINIT (50,2), NINI (2)
      REAL RNDM(4)
      SAVE AMMOLD, PFINIT, TKINIT, TSINIT, TRINIT, NINI
      DATA AMMOLD / 0.9382796D+00, 0.9395731D+00 /
*
      LINIT = .FALSE.
      IF ( NINI (IT) .GT. 0 ) THEN
         PKFRMI = PFINIT (NINI(IT),IT)
         TKIN   = TKINIT (NINI(IT),IT)
         TSTRCK = TSINIT (NINI(IT),IT)
         TKRECL = TRINIT (NINI(IT),IT)
         NINI  (IT) = NINI (IT) - 1
         EXSOP (IT) = 0.D+00
         GO TO 3500
      ELSE
         GO TO 200
      END IF
      ENTRY RBKINI ( IT, LZEROI, EXSOP, TKIN, TSTRCK,
     &               PSTRCK, ARECL, TKRECL )
      LINIT = .TRUE.
      IF ( LZEROI ) THEN
         NINI (1) = 0
         NINI (2) = 0
         RETURN
      END IF
  200 CONTINUE
*  In this version the low energy component has been suppressed, since
*  they are now produced by the evaporation model, A. Ferrari.
*  The parameters needed for the sampling have been already set in
*  Corrin:
      EXSOP (IT) = 0.D+00
*  Sample the Fermi momentum
      ESLPFF = ESLOPE (IT)
      EXUPFF = EXUPNU (IT)
      EKUPFF = EKUPNU (IT)
      EXDWFF = EXMNNU (IT)
      EKDWFF = EKMNNU (IT)
      ERCLFF = ERCLAV (IT)
      IRECNT = 0
*  +-------------------------------------------------------------------*
*  |
 100  CONTINUE
*  |  Sample the Fermi momentum
         CALL GRNDM(RNDM,4)
         PKFRMI = PFRMMX (IT) * MAX ( RNDM (1), RNDM (2),
     &            RNDM (3) )
         PKFRSQ = PKFRMI * PKFRMI
         TKIN   = - ESLPFF * LOG ( EXDWFF - RNDM (4) * (
     &              EXDWFF - EXUPFF ) )
         TKRECL = 0.5D+00 * PKFRSQ / ( AMUC12 * ARECL ) * ( 1.D+00 -
     &            0.25D+00 * PKFRSQ / ( ARECL**2 * AMUCSQ ) )
         TSTRCK = SQRT ( AMNUSQ (IT) + PKFRSQ ) + TKIN - AMNUCL (IT)
     &          - V0WELL (IT) - TKRECL - EBNDNG (IT)
*  |  +----------------------------------------------------------------*
*  |  |
         IF ( TSTRCK .LE. VEFFNU (IT) - V0WELL (IT) ) THEN
            EXSOP (IT) = EXSOP (IT) + TKIN
            IRECNT = IRECNT + 1
*  |  |  +-------------------------------------------------------------*
*  |  |  |
            IF ( IRECNT .GT. 10 ) THEN
*  |  |  |  +----------------------------------------------------------*
*  |  |  |  |
               IF ( TKRECL - ERCLFF .GT. EKUPFF - EKDWFF .AND.
     &              IRECNT .LT. 20 ) THEN
                  ERCLFF = ERCLFF + TKRECL - ERCLFF
                  EKUPFF = EKUPFF + TKRECL - ERCLFF
                  EKDWFF = EKDWFF + TKRECL - ERCLFF
                  AHELP  = EXP ( - ( TKRECL - ERCLFF ) / ESLPFF )
                  EXUPFF = EXUPFF * AHELP
                  EXDWFF = EXDWFF * AHELP
*  |  |  |  |
*  |  |  |  +----------------------------------------------------------*
*  |  |  |  |
               ELSE IF ( IRECNT .GT. 15 ) THEN
                  TKRECL = MAX ( TKRECL, ERCLFF )
                  ERCLFF = ERCLFF + TKRECL
                  EKUPFF = EKUPFF + TKRECL - ERCLFF
                  EKDWFF = EKDWFF + TKRECL - ERCLFF
                  AHELP  = EXP ( - TKRECL / ESLPFF )
                  EXUPFF = EXUPFF * AHELP
                  EXDWFF = EXDWFF * AHELP
               END IF
*  |  |  |  |
*  |  |  |  +----------------------------------------------------------*
            END IF
*  |  |  |
*  |  |  +-------------------------------------------------------------*
            GO TO 100
*  +-<|--<--<--<--<--< go to resampling
         END IF
*  |  |
*  |  +----------------------------------------------------------------*
*  +-------------------------------------------------------------------*
*  |
      IF ( LINIT ) THEN
         NINI (IT) = NINI (IT) + 1
         PFINIT (NINI(IT),IT) = PKFRMI
         TKINIT (NINI(IT),IT) = TKIN
         TSINIT (NINI(IT),IT) = TSTRCK
         TRINIT (NINI(IT),IT) = TKRECL
         RETURN
      END IF
*  |
*  +-------------------------------------------------------------------*
 3500 CONTINUE
*  | Masses have been updated
      PSTRCK = SQRT ( TSTRCK * ( TSTRCK + 2.D+00 * AMNUCL (IT) ) )
*
********************* Sample the angle ********************************
* Polar angle selection
      ADE=0.090D0*(1.D0+0.081D0*ATO1O3)/TKIN
      DEX=EXP(- PIO2 * PIO2 / ADE)
      AN1=(1.D0-DEX)*ADE/2.D0
      AN2=DEX*PIO2
      AN=AN1+AN2
      AN1=AN1/AN
      CALL GRNDM(RNDM,1)
      IF(RNDM(1).GT.AN1)  GO TO 3
    2 CONTINUE
      CALL GRNDM(RNDM,1)
      DE=SQRT(-ADE*LOG(1.D0-RNDM(1)*(1.D0-DEX)))
      IF(DE.GT.PIO2) GO TO 2
C     WRITE(LUNOUT,*)IT,TO,AMSS,SQAMSS,T,P,DE
      COD = COS (DE)
      SID = SIN (DE)
      RETURN
    3 CONTINUE
      CALL GRNDM(RNDM,1)
      COD = - RNDM(1)
      SID = SQRT ( (1.D0 + COD ) * ( 1.D0 - COD ) )
C     WRITE(LUNOUT,*)IT,TO,AMSS,SQAMSS,T,P,DE
      RETURN
      ENTRY RBKMIN (IT)
      NINI(IT) = NINI(IT)-1
      RETURN
      END
