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
*$ CREATE EVDEEX.FOR
*COPY EVDEEX
*
*=== evdeex ===========================================================*
*
      SUBROUTINE EVDEEX ( WEE )
 
#include "geant321/dblprc.inc"
#include "geant321/dimpar.inc"
#include "geant321/iounit.inc"
*
*----------------------------------------------------------------------*
*                                                                      *
*  Evdeex :  created on 5-10-1990 by A. Ferrari & P. Sala, INFN Milan  *
*                                                                      *
*    Last change  on  28-jan-93   by   Alfredo Ferrari, INFN-Milan     *
*                                                                      *
*    This routine provides a simple model for sampling nuclear deexci- *
*    tation gammas following the evaporation step                      *
*                                                                      *
*----------------------------------------------------------------------*
*
#include "geant321/balanc.inc"
#include "geant321/eva0.inc"
#include "geant321/finuc.inc"
#include "geant321/nucdat.inc"
#include "geant321/parevt.inc"
#include "geant321/resnuc.inc"
*
*----------------------------------------------------------------------*
*  Entering the routine we have:                                       *
*            Ammres is the atomic mass of the residual nucleus         *
*            Ibres (Anow) its mass number                              *
*            Icres (Znow) its atomic number                            *
*            Eres         its total energy                             *
*            Ptres        its momentum                                 *
*            Pxres        x-component of the momentum                  *
*            Pyres        y-component of the momentum                  *
*            Pzres        z-component of the momentum                  *
*            Tvrecl       kinetic energy                               *
*            Tvcms        excitation energy                            *
*----------------------------------------------------------------------*
*
      PARAMETER ( C0M1E1 = 0.306 D+00 )
      PARAMETER ( C0E2E1 = 7.1   D-01 )
      PARAMETER ( HNDFE1 = 5.0   D-03 )
      PARAMETER ( HNDFM1 = 5.0   D-02 )
      PARAMETER ( HNDFE2 = 1.0   D+01 )
*
      DIMENSION COSGAM (3)
      REAL RNDM(2)
*
      IDEEXG = 0
      IF ( IBRES .LE. 0 .OR. TVCMS .LE. GAMMIN ) THEN
         TVCMS = 0.D+00
         RETURN
      END IF
      ECHCK  = ERES - AMMRES
      ECHCK0 = ECHCK
      IBHELP = IBRES / 2
      IF ( IBHELP * 2 .LT. IBRES ) THEN
         CALL EEXLVL ( IBRES, ICRES, DELTA, EEX2ND, EEXDUM )
         IPAR  = 1
         JPAR  = 1
      ELSE
         ICHELP = ICRES / 2
         JPAR   = 0
         IF ( ICHELP * 2 .LT. ICRES ) THEN
            CALL EEXLVL ( IBRES, ICRES, DELTA, EEX1ST, EEXDUM )
            IPAR  = 1
         ELSE
            CALL EEXLVL ( IBRES, ICRES, EEX1ST, DELTA, EEXDUM )
            IPAR  = 2
         END IF
      END IF
      DELPAI = MIN ( EEXDUM, DELTA )
      RNUCL  = R0NUCL * RMASS (IBRES)
      AINERM = 0.24D+00 * AMMRES * RNUCL * RNUCL
      ROTEN0 = 0.5D+00 * PLABRC * PLABRC / AINERM
      ENMIN  = MAX ( DELTA, TWOTWO * ROTEN0 )
      RNMASS = AMMRES + TVCMS
      UMO = RNMASS
      GAMCM  = ERES  / RNMASS
      ETAX   = PXRES / RNMASS
      ETAY   = PYRES / RNMASS
      ETAZ   = PZRES / RNMASS
      IF ( TVCMS .LE. ENMIN .OR. IBRES .LE. 4 ) GO TO 3000
      AHELP  = HNDFE1 * RMASS (IBRES) * RMASS (IBRES)
      CFM1E1 = C0M1E1 * HNDFM1 / AHELP
      CFE2E1 = C0E2E1 * HNDFE2 / AHELP
 1000 CONTINUE
         UMEV   = 1.D+03 * TVCMS
         ASMALL = GETA ( UMEV  , ICRES , IBRES-ICRES, ILVMOD, ISDUM,
     &                   AOGMAX, AOGMIN )
         TEMPSQ = 1.D-03 * ( TVCMS - DELPAI ) / ASMALL
         TEMPER = SQRT ( TEMPSQ )
         HHH    = TVCMS / TEMPER
         HHHSQ  = HHH * HHH
         IF(HHH.GT.88) THEN
            AHELP=0.0
         ELSE
            AHELP  = EXP ( - HHH )
         ENDIF
         BRE1M1 = 6.D+00 - AHELP * ( HHHSQ * HHH + 3.D+00 * HHHSQ
     &          + 6.D+00 * HHH )
         BRE2   = 20.D+00 * BRE1M1 - AHELP * ( HHHSQ * HHHSQ * HHH
     &          + 5.D+00 * HHHSQ * HHHSQ )
         BRE1M1 = BRE1M1 * ( 1.D+00 + CFM1E1 )
         BRE2   = BRE2 * TEMPSQ * CFE2E1
         CALL GRNDM(RNDM,1)
         IF ( RNDM (1) .LT. BRE1M1 / ( BRE1M1 + BRE2 ) ) THEN
            LMULT = 1
         ELSE
            LMULT = 2
         END IF
         LEXPN = 2 * LMULT + 1
         DDLEXP = LEXPN
         XDISMX = MIN ( DDLEXP, HHH )
         IF(XDISMX.GT.88) THEN
            DISMX=0.
         ELSE
            DISMX  = XDISMX ** LEXPN * EXP ( - XDISMX )
         ENDIF
 2000    CONTINUE
            CALL GRNDM(RNDM,2)
            XTENT = RNDM (1) * HHH
            IF(XTENT.GT.88) THEN
               FREJE = 0.
            ELSE
               FREJE = XTENT ** LEXPN * EXP ( - XTENT ) / DISMX
            ENDIF
         IF ( RNDM (2) .GE. FREJE ) GO TO 2000
         ENERG0 = XTENT * TEMPER
         TVCMS  = TVCMS - ENERG0
         RNMASS = AMMRES + TVCMS
         CALL RACO ( COSGAM (1), COSGAM (2), COSGAM (3) )
         ERNCM = 0.5D+00 * ( UMO * UMO + RNMASS * RNMASS ) / UMO
         EEGCM = UMO - ERNCM
         PCMS  = EEGCM
         PCMSX = PCMS * COSGAM (1)
         PCMSY = PCMS * COSGAM (2)
         PCMSZ = PCMS * COSGAM (3)
         ETAPCM = PCMSX * ETAX + PCMSY * ETAY + PCMSZ * ETAZ
         ENERG = GAMCM * EEGCM + ETAPCM
         PHELP = ETAPCM / ( GAMCM + 1.D+00 ) + EEGCM
         PLBGX = PCMSX + ETAX * PHELP
         PLBGY = PCMSY + ETAY * PHELP
         PLBGZ = PCMSZ + ETAZ * PHELP
         ERES  = GAMCM * ERNCM - ETAPCM
         EKRES = ERES - RNMASS
         PHELP = - ETAPCM / ( GAMCM + 1.D+00 ) + ERNCM
         PXRES = - PCMSX + ETAX * PHELP
         PYRES = - PCMSY + ETAY * PHELP
         PZRES = - PCMSZ + ETAZ * PHELP
         ECHCK  = ECHCK - ENERG
         IDEEXG = IDEEXG + 1
         NP     = NP + 1
         WEI   (NP) = WEE
         KPART (NP) = 7
         TKI   (NP) = ENERG
         PLR   (NP) = ENERG
         CXR   (NP) = PLBGX / ENERG
         CYR   (NP) = PLBGY / ENERG
         CZR   (NP) = PLBGZ / ENERG
         GAMCM  = ERES  / RNMASS
         ETAX   = PXRES / RNMASS
         ETAY   = PYRES / RNMASS
         ETAZ   = PZRES / RNMASS
         UMO    = RNMASS
      IF ( TVCMS .GT. ENMIN ) GO TO 1000
      IF ( NP .GE. MXP ) THEN
         WRITE ( LUNOUT, * )' **** Finuc overflow in Evdeex,',
     &                      ' program stopped ****'
         WRITE ( LUNERR, * )' **** Finuc overflow in Evdeex,',
     &                      ' program stopped ****'
         STOP
      END IF
 3000 CONTINUE
         IF ( TVCMS .LE. ( 6 + 2 * JPAR ) * ROTEN0 ) THEN
            ENERG0 = TVCMS
            TVCMS  = 0.D+00
            CALL RACO ( COSGAM (1), COSGAM (2), COSGAM (3) )
            ERNCM = 0.5D+00 * ( UMO * UMO + AMMRES * AMMRES ) / UMO
            EEGCM = UMO - ERNCM
            PCMS  = EEGCM
            PCMSX = PCMS * COSGAM (1)
            PCMSY = PCMS * COSGAM (2)
            PCMSZ = PCMS * COSGAM (3)
            ETAPCM = PCMSX * ETAX + PCMSY * ETAY + PCMSZ * ETAZ
            ENERG = GAMCM * EEGCM + ETAPCM
            PHELP = ETAPCM / ( GAMCM + 1.D+00 ) + EEGCM
            PLBGX = PCMSX + ETAX * PHELP
            PLBGY = PCMSY + ETAY * PHELP
            PLBGZ = PCMSZ + ETAZ * PHELP
            ERES  = GAMCM * ERNCM - ETAPCM
            EKRES = ERES - AMMRES
            PHELP = - ETAPCM / ( GAMCM + 1.D+00 ) + ERNCM
            PXRES = - PCMSX + ETAX * PHELP
            PYRES = - PCMSY + ETAY * PHELP
            PZRES = - PCMSZ + ETAZ * PHELP
            ECHCK  = ECHCK - ENERG
            IDEEXG = IDEEXG + 1
            NP     = NP + 1
            WEI   (NP) = WEE
            KPART (NP) = 7
            TKI   (NP) = ENERG
            PLR   (NP) = ENERG
            CXR   (NP) = PLBGX / ENERG
            CYR   (NP) = PLBGY / ENERG
            CZR   (NP) = PLBGZ / ENERG
         ELSE
            AIAMOM = SQRT ( 0.25D+00 + TVCMS / ROTEN0
     &             + 0.75D+00 * JPAR ) - 0.25D+00
            IF ( IPAR .EQ. 1 ) THEN
               IF ( JPAR .EQ. 1 ) THEN
                  JAMIN = 2 * INT ( AIAMOM ) + 1
                  IAMIN = JAMIN / 2
                  IF ( MOD ( IAMIN, 2 ) .GT. 0 ) THEN
                     EGROUN = 3.75D+00 * ROTEN0
                  ELSE
                     EGROUN = 0.75D+00 * ROTEN0
                  END IF
               ELSE
                  IAMIN = INT ( AIAMOM )
                  JAMIN = 2 * IAMIN
                  IF ( MOD ( IAMIN, 2 ) .GT. 0 ) THEN
                     EGROUN = 2.D+00 * ROTEN0
                  ELSE
                     EGROUN = 0.D+00
                  END IF
               END IF
            ELSE
               IAMIN  = INT ( AIAMOM ) / IPAR
               IAMIN  = IAMIN * IPAR
               JAMIN  = 2 * IAMIN
               EGROUN = 0.D+00
            END IF
            DELTAE = TVCMS + EGROUN - 0.25D+00 * ROTEN0 * JAMIN *
     &             ( JAMIN + 2 )
            DO 4000 JAMOM = JAMIN, 4, -4
               ENERG0 = ROTEN0 * ( 2 * JAMOM - 2 )
     &                + DELTAE
               DELTAE = 0.D+00
               TVCMS  = TVCMS - ENERG0
               RNMASS = AMMRES + TVCMS
               CALL RACO ( COSGAM (1), COSGAM (2), COSGAM (3) )
               ERNCM = 0.5D+00 * ( UMO * UMO + RNMASS * RNMASS ) / UMO
               EEGCM = UMO - ERNCM
               PCMS  = EEGCM
               PCMSX = PCMS * COSGAM (1)
               PCMSY = PCMS * COSGAM (2)
               PCMSZ = PCMS * COSGAM (3)
               ETAPCM = PCMSX * ETAX + PCMSY * ETAY + PCMSZ * ETAZ
               ENERG = GAMCM * EEGCM + ETAPCM
               PHELP = ETAPCM / ( GAMCM + 1.D+00 ) + EEGCM
               PLBGX = PCMSX + ETAX * PHELP
               PLBGY = PCMSY + ETAY * PHELP
               PLBGZ = PCMSZ + ETAZ * PHELP
               ERES  = GAMCM * ERNCM - ETAPCM
               EKRES = ERES - RNMASS
               PHELP = - ETAPCM / ( GAMCM + 1.D+00 ) + ERNCM
               PXRES = - PCMSX + ETAX * PHELP
               PYRES = - PCMSY + ETAY * PHELP
               PZRES = - PCMSZ + ETAZ * PHELP
               ECHCK  = ECHCK - ENERG
               IDEEXG = IDEEXG + 1
               NP     = NP + 1
               WEI   (NP) = WEE
               KPART (NP) = 7
               TKI   (NP) = ENERG
               PLR   (NP) = ENERG
               CXR   (NP) = PLBGX / ENERG
               CYR   (NP) = PLBGY / ENERG
               CZR   (NP) = PLBGZ / ENERG
               GAMCM  = ERES / RNMASS
               ETAX   = PXRES / RNMASS
               ETAY   = PYRES / RNMASS
               ETAZ   = PZRES / RNMASS
               UMO    = RNMASS
 4000       CONTINUE
         END IF
      ECHCK  = ECHCK - EKRES
      IF ( ABS ( ECHCK ) .GT. 1.D-7 * ECHCK0 ) THEN
         WRITE ( LUNERR, * )' **** No energy conservation in Evdeex',
     &                      ' ****', ECHCK, IDEEXG
      END IF
      TVCMS = 0.D+00
      IF ( NP .GT. MXP ) THEN
         WRITE ( LUNOUT, * )' **** Finuc overflow in Evdeex,',
     &                      ' program stopped ****'
         WRITE ( LUNERR, * )' **** Finuc overflow in Evdeex,',
     &                      ' program stopped ****'
         STOP
      END IF
      TVRECL = EKRES
      P2RES = PXRES * PXRES + PYRES * PYRES + PZRES * PZRES
      PTRES = SQRT (P2RES)
*=== End of subroutine Evdeex =========================================*
      RETURN
      END
