*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:19:56  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.42  by  S.Giani
*-- Author :
*$ CREATE EVEVAP.FOR
*COPY EVEVAP
*
*=== evevap ===========================================================*
*
      SUBROUTINE EVEVAP ( WEE )
 
#include "geant321/dblprc.inc"
#include "geant321/dimpar.inc"
#include "geant321/iounit.inc"
*
*----------------------------------------------------------------------*
*                                                                      *
*  EVent EVAPoration: this routine is used to steer both the evapora-  *
*  tion, the high energy fission, possibly a future fragmentation      *
*  and the gamma deexcitation routines                                 *
*                                                                      *
*  Created  on  15  may  1991   by   Alfredo Ferrari & Paola Sala      *
*                                             INFN - Milan             *
*                                                                      *
*  Last change  on 19-apr-93    By   Alfredo Ferrari, INFN - Milan     *
*                                                                      *
*----------------------------------------------------------------------*
*
#include "geant321/balanc.inc"
#include "geant321/eva1.inc"
#include "geant321/fheavy.inc"
#include "geant321/finuc.inc"
#include "geant321/hetc5.inc"
#include "geant321/hetc7.inc"
#include "geant321/hettp.inc"
#include "geant321/higfis.inc"
#include "geant321/labcos.inc"
#include "geant321/nucdat.inc"
#include "geant321/parevt.inc"
#include "geant321/part.inc"
#include "geant321/resnuc.inc"
*
      PARAMETER ( AMUMEV = 1.D+03 * AMUAMU )
*
      COMMON /FKEVNT/ LNUCRI, LHADRI
      LOGICAL LNUCRI, LHADRI
*  The initial excitation energy, mass and charge of the nucleus are
*  put into Ex, Apr, Zpr (common Hetc5)
      EX  = MAX ( 1000 * TVCMS, ANGLGB )
      APR = ANOW
      ZPR = ZNOW
*  Reset the fission/fragmentation counter:
      NFISS = 0
*  Ammres is the atomic mass of the residual nucleus
*  Reset accumulators for the energy conservation check (they are only
*  local)
      EOTEST = AMMRES + TVCMS + TVRECL
      ETEVAP = 0.D+00
*  +-------------------------------------------------------------------*
*  |  Set the variables recording the recoil direction of the residual
*  |  nucleus:
      IF ( PTRES .GT. 0.D+00 ) THEN
         COSLBR (1) = PXRES / PTRES
         COSLBR (2) = PYRES / PTRES
         COSLBR (3) = PZRES / PTRES
*  |
*  +-------------------------------------------------------------------*
*  |  It can happen for pion capture for example that ptres=0
*  |  ( it is always 0 if no "direct" particle is emitted )
      ELSE
         COSLBR (1) = 0.D+00
         COSLBR (2) = 0.D+00
         COSLBR (3) = 1.D+00
      END IF
*  |
*  +-------------------------------------------------------------------*
*  The call to getrig is useless, since we actually need no rotation
*     CALL GETRIG ( ZERZER, ZERZER, ONEONE )
      EREC = 1.D+03 * TVRECL
      CALL FKERUP (0)
*  +-------------------------------------------------------------------*
*  |  Check for fission/fragmentation: if it occurred loop back on the
*  |  fission fragments to possibly evaporate further particles:
      IF ( FISINH ) THEN
         LRNFSS = .TRUE.
         FISINH = .FALSE.
         JFISS  = 0
*  |  +----------------------------------------------------------------*
*  |  |  Update the partial counters of evaporated particles
         DO 40 J = 1,6
            NPARTF (J,JFISS) = NPART (J)
            HEVFIS (JFISS)   = HEVSUM
   40    CONTINUE
*  |  |
*  |  +----------------------------------------------------------------*
*  |  +----------------------------------------------------------------*
*  |  |  The following "do" is not structured as a do since Nfiss can
*  |  |  be incremented during evaporation/fragmentation of the
*  |  |  previously generated fragments
   50    CONTINUE
            JFISS  = JFISS + 1
            AMMRES = 1.D-03 * AMFIS (JFISS)
            PTRES  = 1.D-03 * PPFIS (JFISS)
            EREC = EKFIS (JFISS)
            APR  = AFIS  (JFISS)
            ZPR  = ZFIS  (JFISS)
            EX   = MAX ( UFIS (JFISS), ANGLGB )
            COSLBR (1) = COSLFF (1,JFISS)
            COSLBR (2) = COSLFF (2,JFISS)
            COSLBR (3) = COSLFF (3,JFISS)
*  |  |  The call to getrig is useless, since we need no rotation
*           CALL GETRIG ( ZERZER, ZERZER, ONEONE )
            CALL FKERUP (JFISS)
            ANOW = APR
            ZNOW = ZPR
            ICHLP = NINT (ZNOW)
            IBHLP = NINT (ANOW)
*  |  |  +-------------------------------------------------------------*
*  |  |  |  If we enter this branch the present fragment has been
*  |  |  |  completely evaporated without further fragmentation and
*  |  |  |  it is ready for the final gamma deexcitation and for
*  |  |  |  residual nuclei scoring
            IF ( .NOT. FISINH .AND. IBHLP .GT. 0 ) THEN
               AMTFIS (JFISS) = ANOW * AMUMEV +  FKENER ( ANOW, ZNOW )
               UTFIS  (JFISS) = UU
               RECFIS (JFISS) = EREC
               PPTFIS (JFISS) = SQRT ( EREC * ( EREC + TWOTWO
     &                        * ( AMTFIS (JFISS) + UTFIS (JFISS) ) ) )
               ATFIS  (JFISS) = ANOW
               ZTFIS  (JFISS) = ZNOW
               COSLFF (1,JFISS) = COSLBR (1)
               COSLFF (2,JFISS) = COSLBR (2)
               COSLFF (3,JFISS) = COSLBR (3)
               ETEVAP = ETEVAP + 1.D-03 * ( EREC + AMTFIS (JFISS)
     &                + UTFIS (JFISS) )
*  |  |  |
*  |  |  +-------------------------------------------------------------*
*  |  |  |  Fragment furtherly fragmented or completely evaporated into
*  |  |  |  p,n,d,t,3-He and alphas
            ELSE
               FISINH = .FALSE.
               ATFIS  (JFISS) = ZERZER
               ZTFIS  (JFISS) = ZERZER
            END IF
*  |  |  |
*  |  |  +-------------------------------------------------------------*
*  |  |  +-------------------------------------------------------------*
*  |  |  |  Update the partial counters of evaporated particles
            DO 60 J = 1,6
               NPARTF (J,JFISS) = NPART (J)
               HEVFIS (JFISS)   = HEVSUM
   60       CONTINUE
*  |  |  |
*  |  |  +-------------------------------------------------------------*
         IF ( JFISS .LT. NFISS ) GO TO 50
*  |  |
*  |  +----------------------------------------------------------------*
         FISINH = .FALSE.
      END IF
*  |
*  +-------------------------------------------------------------------*
      IEVNEU = NPART (1)
      IEVPRO = NPART (2)
      IEVDEU = NPART (3)
      IEVTRI = NPART (4)
      IEV3HE = NPART (5)
      IEV4HE = NPART (6)
      IEVAPL = IEVNEU + IEVPRO
      IEVAPH = IEVDEU + IEVTRI + IEV3HE + IEV4HE
*  +-------------------------------------------------------------------*
*  |              Add to the secondary stack the evaporated neutrons
      DO 100 IP = 1, NPART (1)
         NP = NP + 1
         KPART (NP) = 8
         TKI   (NP) = 1.D-03 * EPART ( IP, 1 )
         WEI   (NP) = WEE
         CXR   (NP) = COSEVP ( 1, IP, 1 )
         CYR   (NP) = COSEVP ( 2, IP, 1 )
         CZR   (NP) = COSEVP ( 3, IP, 1 )
         PLR   (NP) = SQRT ( TKI (NP) * ( TKI (NP) + 2.D+00 * AM (8) ) )
         ETEVAP = ETEVAP + TKI (NP) + AMHEAV (1)
  100 CONTINUE
*  |
*  +-------------------------------------------------------------------*
 
*  +-------------------------------------------------------------------*
*  |              Add to the secondary stack the evaporated protons
      DO 200 IP = 1, NPART (2)
         NP = NP + 1
         KPART (NP) = 1
         TKI   (NP) = 1.D-03 * EPART ( IP, 2 )
         WEI   (NP) = WEE
         CXR   (NP) = COSEVP ( 1, IP, 2 )
         CYR   (NP) = COSEVP ( 2, IP, 2 )
         CZR   (NP) = COSEVP ( 3, IP, 2 )
         PLR   (NP) = SQRT ( TKI (NP) * ( TKI (NP) + 2.D+00 * AM (1) ) )
         ETEVAP = ETEVAP + TKI (NP) + AMHEAV (2)
  200 CONTINUE
*  |
*  +-------------------------------------------------------------------*
 
*  +-------------------------------------------------------------------*
*  |         Add to the heavy stack the other evaporated (if requested)
      IF ( LHEAVY ) THEN
         NPHEAV = 0
*  |  +----------------------------------------------------------------*
*  |  |  Loop over the particle types:
         DO 400 JP = 3, 6
*  |  |  +-------------------------------------------------------------*
*  |  |  |
            DO 300 IP = 1, NPART (JP)
               NPHEAV = NPHEAV + 1
               KHEAVY (NPHEAV) = JP
               TKHEAV (NPHEAV) = 1.D-03 * EPART ( IP, JP )
               WHEAVY (NPHEAV) = WEE
               CXHEAV (NPHEAV) = COSEVP ( 1, IP, JP )
               CYHEAV (NPHEAV) = COSEVP ( 2, IP, JP )
               CZHEAV (NPHEAV) = COSEVP ( 3, IP, JP )
               PHEAVY (NPHEAV) = SQRT ( ( TKHEAV (NPHEAV) + TWOTWO
     &                         * AMHEAV (JP) ) * TKHEAV (NPHEAV) )
               ETEVAP = ETEVAP + TKHEAV (NPHEAV) + AMHEAV (JP)
  300       CONTINUE
*  |  |  |
*  |  |  +-------------------------------------------------------------*
  400    CONTINUE
*  |  |
*  |  +----------------------------------------------------------------*
*  |
*  +-------------------------------------------------------------------*
*  |
      ELSE
         NPHEAV = 0
         ETEVAP = ETEVAP + 1.D-03 * HEVSUM + IEVDEU * AMHEAV (3)
     &          + IEVTRI * AMHEAV (4)
     &          + IEV3HE * AMHEAV (5)
     &          + IEV4HE * AMHEAV (6)
      END IF
*  |
*  +-------------------------------------------------------------------*
*  +-------------------------------------------------------------------*
*  |  Fission and/or fragmentation occurred:
      IF ( LRNFSS ) THEN
         TVHEAV = 1.D-03 * HEVSUM
         IF ( ABS ( ETEVAP - EOTEST )/ EOTEST .GT. 1.D-07 ) THEN
            WRITE ( LUNOUT, * )
     &            ' Evevap_fis: failure in energy conservation!!',
     &                        ETEVAP, EOTEST
            WRITE ( LUNERR, * )
     &            ' Evevap_fis: failure in energy conservation!!',
     &                        ETEVAP, EOTEST
         END IF
         TVCHLP = ZERZER
         IDEHLP = 0
*  |  +----------------------------------------------------------------*
*  |  |  Loop on fission/fragmentation fragments
         DO 5000 JFISS = 1, NFISS
            ANOW  = ATFIS (JFISS)
            ZNOW  = ZTFIS (JFISS)
            IBRES = NINT ( ANOW )
            ICRES = NINT ( ZNOW )
*  |  |  +-------------------------------------------------------------*
*  |  |  |  Check the residual nucleus:
            IF ( IBRES .EQ. 0 ) THEN
               AMMRES = ZERZER
               TVCMS  = ZERZER
               TVRECL = ZERZER
               PTRES  = ZERZER
               PXRES  = ZERZER
               PYRES  = ZERZER
               PZRES  = ZERZER
               ERES   = ZERZER
*  |  |  |
*  |  |  +-------------------------------------------------------------*
*  |  |  |  real fragment:
            ELSE
               AMMRES = 1.D-03 * AMTFIS (JFISS)
               TVCMS  = 1.D-03 * UTFIS  (JFISS)
               TVRECL = 1.D-03 * RECFIS (JFISS)
               PTRES  = 1.D-03 * PPTFIS (JFISS)
               PXRES  = PTRES * COSLFF (1,JFISS)
               PYRES  = PTRES * COSLFF (2,JFISS)
               PZRES  = PTRES * COSLFF (3,JFISS)
               ERES   = AMMRES + TVCMS + TVRECL
               EKRES  = TVRECL
            END IF
*  |  |  |
*  |  |  +-------------------------------------------------------------*
*  |  |  +-------------------------------------------------------------*
*  |  |  |  Check if the deexcitation module have to be called
            IF ( LDEEXG ) THEN
               IDEEXG = 0
               CALL EVDEEX ( WEE )
               IDEHLP = IDEHLP + IDEEXG
*  |  |  |
*  |  |  +-------------------------------------------------------------*
*  |  |  |
            ELSE
               TVCHLP = TVCHLP + TVCMS
            END IF
*  |  |  |
*  |  |  +-------------------------------------------------------------*
*  |  |  +-------------------------------------------------------------*
*  |  |  |  Check if fission fragments have to be put on stack
            IF ( LHEAVY .AND. IBRES .GT. 0 ) THEN
               NPHEAV = NPHEAV + 1
               TKHEAV (NPHEAV) = EKRES
               PHEAVY (NPHEAV) = PTRES
               CXHEAV (NPHEAV) = PXRES / PTRES
               CYHEAV (NPHEAV) = PYRES / PTRES
               CZHEAV (NPHEAV) = PZRES / PTRES
               WHEAVY (NPHEAV) = WEE
               KHEAVY (NPHEAV) = 6 + JFISS
               AMHEAV (KHEAVY(NPHEAV)) = AMMRES
               IBHEAV (KHEAVY(NPHEAV)) = IBRES
               ICHEAV (KHEAVY(NPHEAV)) = ICRES
            END IF
*  |  |  |
*  |  |  +-------------------------------------------------------------*
            TVHEAV = TVHEAV + TVRECL
 5000    CONTINUE
*  |  |
*  |  +----------------------------------------------------------------*
         IDEEXG = IDEHLP
         TVCMS  = TVCHLP
         ANOW   = ZERZER
         ZNOW   = ZERZER
         IBRES  = 0
         ICRES  = 0
         AMMRES = ZERZER
         TVRECL = ZERZER
         PTRES  = ZERZER
         PXRES  = ZERZER
         PYRES  = ZERZER
         PZRES  = ZERZER
         ERES   = ZERZER
*  |
*  +-------------------------------------------------------------------*
*  |  Normal evaporation:
      ELSE
         ANOW  = APR
         ZNOW  = ZPR
         IBRES = NINT ( ANOW )
         ICRES = NINT ( ZNOW )
*  |  Ammres is the atomic mass of the residual nucleus
*  |  +----------------------------------------------------------------*
*  |  |  Check the residual nucleus:
         IF ( IBRES .EQ. 0 ) THEN
            AMMRES = ZERZER
            TVCMS  = ZERZER
            TVRECL = ZERZER
            PTRES  = ZERZER
            PXRES  = ZERZER
            PYRES  = ZERZER
            PZRES  = ZERZER
            ERES   = ZERZER
*  |  |
*  |  +----------------------------------------------------------------*
*  |  |
         ELSE
            AMMRES = ANOW * AMUAMU + 1.D-03 * FKENER ( ANOW, ZNOW )
            TVCMS  = 1.D-03 * UU
            TVRECL = 1.D-03 * EREC
            PTRES  = SQRT ( TVRECL * ( TVRECL + 2.D+00 * ( AMMRES +
     &                      TVCMS ) ) )
            PXRES  = PTRES * COSLBR (1)
            PYRES  = PTRES * COSLBR (2)
            PZRES  = PTRES * COSLBR (3)
            ERES   = AMMRES + TVCMS + TVRECL
            EKRES  = TVRECL
         END IF
*  |  |
*  |  +----------------------------------------------------------------*
         TVHEAV = 1.D-03 * HEVSUM
         ETEVAP = ETEVAP + ERES
         IF ( ABS ( ETEVAP - EOTEST )/ EOTEST .GT. 1.D-07 ) THEN
            WRITE ( LUNOUT, * )
     &            ' Evevap: failure in energy conservation!!',
     &                        ETEVAP, EOTEST
            WRITE ( LUNERR, * )
     &            ' Evevap: failure in energy conservation!!',
     &                        ETEVAP, EOTEST
         END IF
*  |   Check if the deexcitation module have to be called
         IF ( LDEEXG ) CALL EVDEEX ( WEE )
      END IF
*  |
*  +-------------------------------------------------------------------*
      RETURN
*=== End of subroutine Evevap =========================================*
      END
