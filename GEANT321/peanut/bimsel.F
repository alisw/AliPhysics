*
* $Id$
*
* $Log$
* Revision 1.2  1997/10/17 10:00:03  mclareni
* Remove SAVE statement for NT
*
* Revision 1.1.1.1  1995/10/24 10:22:00  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.45  by  S.Giani
*-- Author :
*$ CREATE BIMSEL.FOR
*COPY BIMSEL
*
*=== bimsel ===========================================================*
*
      SUBROUTINE BIMSEL ( JPROJ, TXX, TYY, TZZ, LBCHCK )
 
#include "geant321/dblprc.inc"
#include "geant321/dimpar.inc"
#include "geant321/iounit.inc"
*
*----------------------------------------------------------------------*
*----------------------------------------------------------------------*
*
#include "geant321/balanc.inc"
#include "geant321/eva0.inc"
#include "geant321/nucdat.inc"
#include "geant321/nucgeo.inc"
#include "geant321/part.inc"
#include "geant321/parevt.inc"
#include "geant321/resnuc.inc"
*
      PARAMETER ( FEFFEC = 1.518066780142162 D+00 )
      PARAMETER ( BETMAX = 0.4 D+00 )
*
      REAL RNDM(2)
*
      LOGICAL LBCHCK, LFERMI, LLMDBR
*
      SAVE ABTAR , ZZTAR , SIGMP0, SIGMN0,
     &     AMNHLP, RHOBIM, RPRONU, RADPRP, RADPRN, DSKRED, RHRUSF,
     &     AUSFL , ZUSFL , BNDSAV, RADHLP, BFCHLP, BIMCLM, PRCOLP,
     &     PRCOLN, IBTOLD, ICTOLD, KPROJ , NTRIAL, ITFRMI
      SAVE SBHAL0, SBSKI0, SBCEN0, SBCEN1, SBSKI1, SBHAL1
#ifndef CERNLIB_WINNT
      SAVE
#endif
      DATA IBTOLD, ICTOLD / 2*0 /
*
      KPROJ  = JPROJ
      AUSFL  = IBTAR
      ZUSFL  = ICHTAR
      RHRUSF = 1.D+00
      BEPROJ = PNUCCO / ( EKECON + AM (KPROJ) )
      CXIMPC = TXX
      CYIMPC = TYY
      CZIMPC = TZZ
      NTRIAL = 0
      RHOBIM = - AINFNT
      IF ( KPROJ .EQ. 1 .OR. KPROJ .EQ. 8 ) THEN
         IPWELL = 1 + KPROJ / 8
         WLLRED = 1.D+00
         BNDNUC = BNENRG (IPWELL)
      ELSE
         IPWELL = 0
         IF ( IBAR (KPROJ) .EQ. 0 ) THEN
            IF ( KPROJ .LE. 11 ) THEN
               WLLRED = 0.D+00
               BNDNUC = 0.D+00
            ELSE
               WLLRED = POTMES
               BNDNUC = BNENRG (3)
            END IF
         ELSE
            WLLRED = POTBAR
            BNDNUC = BNENRG (3)
         END IF
      END IF
      BNDSAV = BNDNUC
      IF ( IBAR (KPROJ) .NE. 0 ) THEN
         RPRONU = 1.D+00
      ELSE IF ( KPROJ .NE. 7 ) THEN
         RPRONU = 0.8164965809277260D+00
      ELSE
         RPRONU = 0.D+00
      END IF
      IF ( LBCHCK ) THEN
         LFERMI = .FALSE.
         EKESIG = EKECON
         PPRSIG = PNUCCO
         CALL SIGFER ( KPTOIP (KPROJ), EKESIG, PPRSIG, LFERMI )
         PRCOLP = ZUSFL / AUSFL * SIGMAP
         PRCOLN = ( AUSFL - ZUSFL ) / AUSFL * SIGMAN
         SIGMAA = PRCOLP + PRCOLN
         PRCOLP = PRCOLP / SIGMAA
         PRCOLN = 1.D+00 - PRCOLP
      END IF
      IF ( RPRONU .GT. ANGLGB ) THEN
         IF ( LPARWV ) THEN
            LLMDBR = .TRUE.
            TMP102 = 1.D-02
            PDEBRO = MAX ( PNUCCO, TMP102 )
            ALMBAR = PLABRC / PDEBRO
            DEBRLM = 0.5D+00 * ALMBAR
            RADCOR = SQRT ( (RPRONU * RMSPRO)**2 + ALMBAR**2 )
     &             / ( RMSPRO * RPRONU )
         ELSE
            PDEBRO = ( EKECON + BNDNUC ) * ( EKECON + BNDNUC + 2.D+00
     &             * AM (KPROJ) )
            LLMDBR = .FALSE.
            DEBRLM = 0.D+00
            ALMBAR = 0.D+00
            LLLMAX = -1
            RADCOR = SQRT ( (RPRONU * RMSPRO)**2 + PLABRC**2 / PDEBRO
     &             ) / ( RMSPRO * RPRONU )
         END IF
      ELSE
         RADCOR = 0.D+00
         LLMDBR = .FALSE.
         DEBRLM = 0.D+00
      END IF
      RADCO2 = RADCOR
      RADPRO = MIN ( TWOTWO * RMSPRO * RPRONU * RADCOR, SKGT16 )
      RADPRP = RADPRO
      RADPRN = RADPRO
      IF ( IBTAR .NE. IBTOLD .OR. ICHTAR .NE. ICTOLD ) THEN
         IBTOLD = IBTAR
         ICTOLD = ICHTAR
         ABTAR  = IBTAR
         ZZTAR  = ICHTAR
         AR1O3  = RMASS (IBTAR)
         AMNHLP = 0.5D+00 * ( AMNUCL (1) + AMNUCL (2) )
         HKAP   = ABTAR**2 / ( ZZTAR**2 + ( ABTAR - ZZTAR )**2 )
         HHLP (1) = ( HKAP * ZZTAR )**0.3333333333333333D+00 / AR1O3
         HHLP (2) = ( HKAP * ( ABTAR - ZZTAR ) )
     &            **0.3333333333333333D+00 / AR1O3
         RHOCEN = RHOTAB (IBTAR)
         ALPHAL = ALPTAB (IBTAR)
         RADIU0 = RADTAB (IBTAR)
         SKINDP = SKITAB (IBTAR)
         HALODP = HALTAB (IBTAR)
         RADIU1 = RADIU0 + SKINDP
         RADTOT = RADIU1 + HALODP
         RHOCOR = ONEMNS * RHOCEN
         RHOSKN = ALPHAL * RHOCEN
         PFRCEN (1) = HHLP (1) * PFRTAB (IBTAR)
         PFRCEN (2) = HHLP (2) * PFRTAB (IBTAR)
         RHOAVE = RHATAB (IBTAR)
         PFRAVE = PFATAB (IBTAR)
         EKFAVE = EKATAB (IBTAR)
         OMALHL = 1.D+00 - ALPHAL
         RAD1O2 = RADIU0 + 0.5D+00 * SKINDP / OMALHL
         SKNEFF = SKINDP / OMALHL
         RADSKN = RADIU0 + SKNEFF
         EKFCEN (1) = SQRT ( AMNUSQ (1) + PFRCEN (1)**2 ) - AMNUCL (1)
         EKFCEN (2) = SQRT ( AMNUSQ (2) + PFRCEN (2)**2 ) - AMNUCL (2)
         IF ( PFRCEN (1) .GT. PFRCEN (2) ) THEN
            ITNCMX = 1
         ELSE
            ITNCMX = 2
         END IF
         CALL NCLVST ( IBTAR, ICHTAR )
      END IF
      IF ( IPWELL .LE. 0 ) IPWELL = ITNCMX
      CALL NCLVIN
      IF ( EKECON .LT. 2.D+00 * GAMMIN ) THEN
         EKECON = 0.D+00
         PNUCCO = 0.D+00
         LABRST = .TRUE.
         RADPRO = 0.D+00
         RADSIG = ( RADTOT + DEBRLM ) * BFCLMB
         RADMAX = RADTOT
         LLLMAX = -1
         OPACTY = 2.D+00
         CALL RSTSEL (KPROJ)
         RETURN
      END IF
      RADMAX = RADTOT + RADPRO
      BIMCLM = RDCLMB * BFCLMB
      IF ( LLMDBR ) THEN
         RADHLP = RADMAX
         IF ( RADHLP .LE. RDCLMB ) THEN
            BFCMAX = BFCLMB
            BFCHLP = 0.5D+00 * CLMBBR * RDCLMB / EKECON
         ELSE
            BFCHLP = 0.5D+00 * CLMBBR * RDCLMB / EKECON
            BFCMAX = SQRT ( 1.D+00 - CLMBBR * RDCLMB / EKECON / RADHLP )
         END IF
         BIMMAX = RADHLP * BFCMAX
         LLLMAX = INT ( BIMMAX / ALMBAR )
         RADSIG = ALMBAR * ( LLLMAX + 1.D+00 )
         SIGGEO = PI * RADSIG * RADSIG
      ELSE
         RADHLP = RADTOT + RADPRO + DEBRLM
         IF ( RADHLP .LE. RDCLMB ) THEN
            BFCMAX = BFCLMB
         ELSE
            BFCHLP = 0.5D+00 * CLMBBR * RDCLMB / EKECON
            BFCMAX = SQRT ( 1.D+00 - CLMBBR * RDCLMB / EKECON / RADHLP )
         END IF
         RADSIG = RADHLP * BFCMAX
      END IF
      R0TRAJ = - RADTOT
      R1TRAJ = - R0TRAJ
 4200  CONTINUE
          CALL GRNDM(RNDM,2)
          RPHI1 = 2.D+00 * RNDM (1) - 1.D+00
          RPHI2 = 2.D+00 * RNDM (2) - 1.D+00
          RPHI12 = RPHI1 * RPHI1
          RPHI22 = RPHI2 * RPHI2
          RSQ = RPHI12 + RPHI22
       IF ( RSQ .GT. 1.D+00 ) GO TO 4200
      SINPHI = 2.D+00 * RPHI1 * RPHI2 / RSQ
      COSPHI = ( RPHI12 - RPHI22 ) / RSQ
      SINT02 = 1.D+00 - CZIMPC * CZIMPC
      IF ( SINT02 .LT. ANGLSQ ) THEN
         UBIMPC = COSPHI
         VBIMPC = SINPHI
         WBIMPC = 0.D+00
      ELSE
         SINTH0 = SQRT ( SINT02 )
         SINPH0 = CYIMPC / SINTH0
         COSPH0 = CXIMPC / SINTH0
         UBIMPC = COSPHI * COSPH0 * CZIMPC - SINPHI * SINPH0
         VBIMPC = COSPHI * SINPH0 * CZIMPC + SINPHI * COSPH0
         WBIMPC = - COSPHI * SINTH0
      END IF
      GO TO 4500
      ENTRY BIMNXT ( LBCHCK )
         IF ( EKECON .LT. 2.D+00 * GAMMIN ) THEN
            LABRST = .TRUE.
            CALL RSTNXT
            RETURN
         END IF
 4300 CONTINUE
         BNDNUC = BNDSAV
         SIGMAP = SIGMP0
         SIGMAN = SIGMN0
 4400    CONTINUE
         CALL GRNDM(RNDM,1)
         ANMFP = - LOG ( 1.D+00 - RNDM (1) ) / DSKRED
         IF ( SBRES * SIGMAA .GT. ANMFP ) THEN
            GO TO 6000
         END IF
 4500 CONTINUE
 5000 CONTINUE
         SBUSED = 0.D+00
         NTRIAL = NTRIAL + 1
         IF ( LLMDBR ) THEN
            ALLMAX = LLLMAX + 1.D+00
            CALL GRNDM(RNDM,2)
            RNDLLL = ALLMAX * MAX ( RNDM (1), RNDM (2) )
            LLLACT = INT (RNDLLL)
            BIMPTR = RNDLLL * ALMBAR
            BIMPTR = ABS (BIMPTR)
            IF ( BIMPTR .LE. BIMCLM ) THEN
               BFCEFF = BFCLMB
            ELSE
               HLPHLP = BFCHLP / BIMPTR
               BFCEFF = 1.D+00 / ( HLPHLP + SQRT ( HLPHLP * HLPHLP
     &                + 1.D+00 ) )
            END IF
            BIMPTR = BIMPTR / BFCEFF
            IF ( BIMPTR .GT. RADHLP ) GO TO 5000
         ELSE
            CALL GRNDM(RNDM,2)
            BIMPTR = RADSIG * MAX ( RNDM (1), RNDM (2) )
            IF ( BIMPTR .LE. BIMCLM ) THEN
               BFCEFF = BFCLMB
            ELSE
               HLPHLP = BFCHLP / BIMPTR
               BFCEFF = 1.D+00 / ( HLPHLP + SQRT ( HLPHLP * HLPHLP
     &                + 1.D+00 ) )
            END IF
            BIMPTR = BIMPTR / BFCEFF
         END IF
         BIMMEM = BIMPTR
         IF ( BIMPTR .GT. RADTOT - RADPRO ) THEN
            BIMPCT = 0.5D+00 * ( RADTOT + BIMPTR - RADPRO )
            IF ( BIMPTR .GE. RADTOT ) THEN
               X1 = BIMPTR - RADTOT
               ANGRED = ACOS ( 2.D+00 * X1 / ( RADPRO + X1 ) ) / PI
               X1 = X1 / ( R0PROT * RPRONU * RADCO2 )
               DSKRED = ( 0.5D+00 * X1 * X1 + X1 + 1.D+00 ) * EXP (-X1)
     &                * ANGRED
            ELSE
               X1 = RADPRO + BIMPTR - RADTOT
               ANGRED = ACOS ( 2.D+00 * X1 / ( RADPRO + X1 ) ) / PI
               X1 = X1 / ( R0PROT * RPRONU * RADCO2 )
               DSKRED = 1.D+00 - ( 0.5D+00 * X1 * X1 + X1 + 1.D+00 )
     &                * EXP (-X1) * ANGRED
            END IF
         ELSE
            DSKRED = 1.D+00
            BIMPCT = BIMPTR
         END IF
         IF ( .NOT. LBCHCK ) THEN
            RHOSAV = RHOBIM
            RHOBIM = FRHONC ( BIMPCT )
            IF ( RHOBIM .EQ. RHOSAV ) GO TO 5500
            PFRBIM = FPFRNC ( RHOBIM, ITNCMX )
            EKFBIM = FEKFNC ( PFRBIM, ITNCMX )
            RHOHLP = FRHONC ( BIMPTR )
            PFRHLP = FPFRNC ( RHOHLP, IPWELL )
            PFRHLP = 0.5D+00 * PFRHLP * PFRHLP / AMNUSQ (IPWELL)
            IF ( BIMPTR .GT. RADTOT ) BNDNUC = BNDNUC * ( 1.D+00
     &         - ( BIMPTR - RADTOT ) / ( RADHLP - RADTOT ) )
            VPRBIM = WLLRED * ( AMNUCL (IPWELL) * PFRHLP
     &             * ( 1.D+00 - 0.5D+00 * PFRHLP ) + BNDNUC )
            LFERMI = .TRUE.
            EKESIG = EKECON
            PPRSIG = PNUCCO
            CALL SIGFER ( KPTOIP (KPROJ), EKESIG, PPRSIG, LFERMI )
            PRCOLP = ZUSFL / AUSFL * SIGMAP
            PRCOLN = ( AUSFL - ZUSFL ) / AUSFL * SIGMAN
            SIGMAA = PRCOLP + PRCOLN
            PRCOLP = PRCOLP / SIGMAA
            PRCOLN = 1.D+00 - PRCOLP
 5500       CONTINUE
         ELSE
            RHOBIM = FRHONC ( BIMPCT )
         END IF
         XBIMPC = UBIMPC * BIMPCT
         YBIMPC = VBIMPC * BIMPCT
         ZBIMPC = WBIMPC * BIMPCT
         CALL GRNDM(RNDM,1)
         ANMFP  = - LOG ( 1.D+00 - RNDM (1) ) / DSKRED
         IF ( BIMPCT .GT. RAD1O2 ) THEN
            SBTTSQ = 4.D+00 * ( RADTOT**2 - BIMPCT**2 ) * RHOBIM**2
            IF ( SBTTSQ .LE. ( ANMFP / SIGMAA )**2 ) GO TO 5000
         END IF
         CALL SBCOMP ( SBHAL0, SBSKI0, SBCEN0, SBCEN1, SBSKI1, SBHAL1 )
         SBTOT  = SBHAL0 + SBSKI0 + SBCEN0 + SBCEN1 + SBSKI1 + SBHAL1
         SBTOT  = RHRUSF * SBTOT
      IF ( SBTOT * SIGMAA .LE. ANMFP ) GO TO 5000
 6000 CONTINUE
      SBUSED = SBUSED * RHRUSF + ANMFP / SIGMAA
      SBRES  = SBTOT  - SBUSED
      SBUSED = SBUSED / RHRUSF
      LELSTC = .TRUE.
      NTARGT = 1
      CALL GRNDM(RNDM,1)
      IF ( RNDM (1) .LT. PRCOLP ) THEN
         KNUCIM = 1
         ITFRMI = 1
      ELSE
         KNUCIM = 8
         ITFRMI = 2
      END IF
      IPRTYP = KPROJ * 10 + KNUCIM
      CALL RSCOMP ( SBHAL0, SBSKI0, SBCEN0, SBCEN1, SBSKI1, SBHAL1 )
      BNDNUC = BNDSAV
      RHOIMP = FRHONC ( ABS (RIMPCT) )
      PFRIMP = FPFRNC ( RHOIMP, ITFRMI )
      EKFIMP = FEKFNC ( PFRIMP, ITFRMI )
      RIMHLP = ABS (RIMPTR)
      RHOIMT = FRHONC ( RIMHLP )
      PFRPRO = FPFRNC ( RHOIMT, IPWELL )
      EKFPRO = FEKFNC ( PFRPRO, IPWELL )
      IF ( RIMHLP .GT. RADTOT ) BNDNUC = BNDNUC * ( 1.D+00 - ( RIMHLP
     &                                 - RADTOT ) / ( RADHLP - RADTOT ))
      VPRWLL = WLLRED * ( EKFPRO + BNDNUC )
      EKEWLL = EKECON + VPRWLL
      EPSWLL = EKEWLL + AM (KPROJ)
      IF ( .NOT. LBCHCK ) THEN
         PPRWLL = SQRT ( EKEWLL * ( EKEWLL + 2.D+00 * AM (KPROJ) ) )
         CALL PHDWLL ( UBIMPC, VBIMPC, WBIMPC )
         PNFRMI = PFNCLV ( ITFRMI, .TRUE. )
         IF ( PNFRMI .LT. -100.D+00 ) GO TO 4400
         CALL RACO ( PXFERM, PYFERM, PZFERM )
         PXFERM = PXFERM * PNFRMI
         PYFERM = PYFERM * PNFRMI
         PZFERM = PZFERM * PNFRMI
         ERES   = EKEWLL + AM (KPROJ) + AM (KNUCIM) + EKFERM
         PXRES  = PXPROJ + PXFERM
         PYRES  = PYPROJ + PYFERM
         PZRES  = PZPROJ + PZFERM
         PTRES2 = PXRES**2 + PYRES**2 + PZRES**2
         UMO2   = ERES**2 - PTRES2
         EKESIG = 0.5D+00 * ( UMO2 - AM (KPROJ)**2 - AM (KNUCIM)**2 )
     &          / AM (KNUCIM) - AM (KPROJ)
         EKFIMP = MAX ( EKFERM, EKFIMP )
      ELSE
         EKESIG = EPSWLL - AM (KPROJ)
      END IF
      PPRSIG = SQRT ( EKESIG * ( EKESIG + 2.D+00 * AM (KPROJ) ) )
      SIGMN0 = SIGMAN
      SIGMP0 = SIGMAP
      LFERMI = .FALSE.
      CALL SIGFER ( KPTOIP (KPROJ), EKESIG, PPRSIG, LFERMI )
      IF ( KNUCIM .EQ. 1 ) THEN
         SIGMAR = SIGMAP / SIGMP0
      ELSE
         SIGMAR = SIGMAN / SIGMN0
      END IF
      SIGMAR = MIN ( SIGMAR, ONEONE )
      CALL GRNDM(RNDM,1)
      RNDREJ = RNDM(1)
      IF ( RNDREJ .GE. SIGMAR ) GO TO 4300
      IF ( LBCHCK ) THEN
         ZITA   = 0.5D+00 * ( EKFIMP + EKFPRO ) / EKEWLL
         IF ( ZITA .LE. 0.5D+00 ) THEN
            PZITA = 1.D+00 - 1.4D+00 * ZITA
         ELSE
            PZITA = 1.D+00 - 1.4D+00 * ZITA + 0.4D+00 * ZITA
     &            * ( 2.D+00 - 1.D+00 / ZITA )**2.5D+00
         END IF
         RNDREJ = RNDREJ / SIGMAR
         IF ( RNDREJ .GE. PZITA ) GO TO 4300
      ELSE
         PZITA = 1.D+00
      END IF
      OPACTY = 1.D+00 / NTRIAL
      RETURN
*=== End of subroutine Bimsel =========================================*
      END
