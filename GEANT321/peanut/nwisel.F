*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:22:01  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.46  by  S.Giani
*-- Author :
*$ CREATE NWISEL.FOR
*COPY NWISEL
*
*=== nwisel ===========================================================*
*
      SUBROUTINE NWISEL ( IKN, LBCHCK )
 
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
#include "geant321/parevt.inc"
#include "geant321/parnuc.inc"
#include "geant321/part.inc"
#include "geant321/resnuc.inc"
*
      PARAMETER ( FEFFEC = 1.518066780142162 D+00 )
      PARAMETER ( BETMAX = 0.4 D+00 )
*
      REAL RNDM(1), RNGAUS, DUMNOR
*
      LOGICAL LBCHCK, LFERMI, LFRMCK, LLMDBR
*
      SAVE XBIMTR, YBIMTR, ZBIMTR, COSTHE, SINTHE, SIGMP0, EKENWI,
     &     SIGMN0, RHOBIM, RPRONU, RADPRP, RADPRN, DSKRED, RHRUSF,
     &     AUSFL , ZUSFL , PRCOLP, PRCOLN, KPROJ , IKPMX , LFRMCK
      SAVE SBHAL0, SBSKI0, SBCEN0, SBCEN1, SBSKI1, SBHAL1
*
      LFRMCK = IKN .LT. 0
      IKPMX  = ABS (IKN)
      KPROJ  = KPNUCL (IKPMX)
      EKENWI = EKECON
      AUSFL  = IBTAR
      ZUSFL  = ICHTAR
      RHRUSF = 1.D+00
      BEPROJ = PNUCCO / ( EKENWI + AM (KPROJ) )
      RHOBIM = - AINFNT
*
      IF ( KPROJ .EQ. 1 .OR. KPROJ .EQ. 8 ) THEN
         IPWELL = 1 + KPROJ / 8
         WLLRED = 1.D+00
      ELSE
         IPWELL = ITNCMX
         IF ( IBAR (KPROJ) .EQ. 0 ) THEN
            IF ( KPROJ .LE. 11 ) THEN
               WLLRED = 0.D+00
            ELSE
               WLLRED = POTMES
            END IF
         ELSE
            WLLRED = POTBAR
         END IF
      END IF
      IF ( IBAR (KPROJ) .NE. 0 ) THEN
         RPRONU = 1.D+00
      ELSE IF ( KPROJ .NE. 7 ) THEN
         RPRONU = 0.8164965809277260D+00
      ELSE
         RPRONU = 0.D+00
      END IF
      IF ( LBCHCK ) THEN
         LFERMI = .FALSE.
         EKESIG = EKENWI
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
            PDEBRO = MAX ( PNUCL (IKPMX), TMP102 )
            DEBRLM = 0.5D+00 * PLABRC / PDEBRO
            RADCOR = 1.D+00
         ELSE
            LLMDBR = .FALSE.
            DEBRLM = 0.D+00
            RADCOR = 1.D+00
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
      CXIMPC = PXNUCL (IKPMX) / PNUCL (IKPMX)
      CYIMPC = PYNUCL (IKPMX) / PNUCL (IKPMX)
      CZIMPC = PZNUCL (IKPMX) / PNUCL (IKPMX)
      COSTHE = ( CXIMPC * XSTNUC (IKPMX) + CYIMPC * YSTNUC (IKPMX)
     &         + CZIMPC * ZSTNUC (IKPMX) ) / RSTNUC (IKPMX)
      SINTHE = SQRT ( 1.D+00 - MIN ( COSTHE * COSTHE, ONEONE ) )
      BIMPTR = RSTNUC (IKPMX) * SINTHE
      IF ( LLMDBR ) THEN
         BIMOLD = BIMPTR
         CALL GRNDM(RNDM,1)
         BIMPTR = BIMPTR + ( 2.D+00 * RNDM (1) - 1.D+00 )
     &          * DEBRLM
         CALL GRANOR ( RNGAUS, DUMNOR)
         HEISEX = DEBRLM * RNGAUS
         DSTMIN = - RSTNUC (IKPMX) * COSTHE
         XBIMTR = XSTNUC (IKPMX) + DSTMIN * CXIMPC
         YBIMTR = YSTNUC (IKPMX) + DSTMIN * CYIMPC
         ZBIMTR = ZSTNUC (IKPMX) + DSTMIN * CZIMPC
         XPARAL = XSTNUC (IKPMX) - XBIMTR
         YPARAL = YSTNUC (IKPMX) - YBIMTR
         ZPARAL = ZSTNUC (IKPMX) - ZBIMTR
         IF ( BIMOLD .GT. ANGLGB ) THEN
            BIMOLD = BIMPTR / BIMOLD
            XBIMTR = XBIMTR * BIMOLD
            YBIMTR = YBIMTR * BIMOLD
            ZBIMTR = ZBIMTR * BIMOLD
         ELSE
            CALL RACO ( CXHLP, CYHLP, CZHLP )
            XBIMTR = BIMPTR * ( CYHLP * CZIMPC - CZHLP * CYIMPC )
            YBIMTR = BIMPTR * ( CZHLP * CXIMPC - CXHLP * CZIMPC )
            ZBIMTR = BIMPTR * ( CXHLP * CYIMPC - CYHLP * CXIMPC )
         END IF
         XSTNUC (IKPMX) = XPARAL + XBIMTR + HEISEX * CXIMPC
         YSTNUC (IKPMX) = YPARAL + YBIMTR + HEISEX * CYIMPC
         ZSTNUC (IKPMX) = ZPARAL + ZBIMTR + HEISEX * CZIMPC
         RSTNUC (IKPMX) = SQRT ( XSTNUC (IKPMX)**2 + YSTNUC (IKPMX)**2
     &                  + ZSTNUC (IKPMX)**2 )
         COSTHE = ( CXIMPC * XSTNUC (IKPMX) + CYIMPC * YSTNUC (IKPMX)
     &          + CZIMPC * ZSTNUC (IKPMX) ) / RSTNUC (IKPMX)
         BIMPTR = ABS (BIMPTR)
         SINTHE = BIMPTR / MAX ( RSTNUC (IKPMX), ANGLGB )
      ELSE
         DSTMIN = - RSTNUC (IKPMX) * COSTHE
         XBIMTR = XSTNUC (IKPMX) + DSTMIN * CXIMPC
         YBIMTR = YSTNUC (IKPMX) + DSTMIN * CYIMPC
         ZBIMTR = ZSTNUC (IKPMX) + DSTMIN * CZIMPC
      END IF
      IF ( COSTHE .GE. 0.D+00 .AND. RSTNUC (IKPMX) .GE. RADTOT
     &   + RADPRO ) GO TO 4500
*
      SBUSED = 0.D+00
      IF ( BIMPTR .GT. RADTOT - RADPRO ) THEN
         BIMPCT = 0.5D+00 * ( RADTOT + BIMPTR - RADPRO )
         TRUFAC = BIMPCT / BIMPTR
         IF ( BIMPTR .GE. RADTOT ) THEN
            X1 = BIMPTR - RADTOT
            IF ( X1 .LE. RADPRO ) THEN
               ANGRED = ACOS ( 2.D+00 * X1 / ( RADPRO + X1 ) ) / PI
               X1 = X1 / ( R0PROT * RPRONU * RADCO2 )
               DSKRED = ( 0.5D+00 * X1 * X1 + X1 + 1.D+00 )
     &                * EXP (-X1) * ANGRED
            END IF
         ELSE
            X1 = RADPRO + BIMPTR - RADTOT
            ANGRED = ACOS ( 2.D+00 * X1 / ( RADPRO + X1 ) ) / PI
            X1 = X1 / ( R0PROT * RPRONU * RADCO2 )
            DSKRED = 1.D+00 - ( 0.5D+00 * X1 * X1 + X1 + 1.D+00 )
     &             * EXP (-X1) * ANGRED
         END IF
      ELSE
         BIMPCT = BIMPTR
         TRUFAC = 1.D+00
         DSKRED = 1.D+00
      END IF
      IF ( BIMPCT .LT. RADTOT ) THEN
         IF ( .NOT. LBCHCK ) THEN
            RHOSAV = RHOBIM
            RHOBIM = FRHONC ( BIMPCT )
            IF ( RHOBIM .EQ. RHOSAV ) GO TO 5500
            PFRBIM = FPFRNC ( RHOBIM, ITNCMX )
            EKFBIM = FEKFNC ( PFRBIM, ITNCMX )
            RHOHLP = FRHONC ( BIMPTR )
            PFRHLP = FPFRNC ( RHOHLP, IPWELL )
            PFRHLP = 0.5D+00 * PFRHLP * PFRHLP / AMNUSQ (IPWELL)
            VPRBIM = WLLRED * ( AMNUCL (IPWELL) * PFRHLP
     &             * ( 1.D+00 - 0.5D+00 * PFRHLP ) + BNDNUC )
            LFERMI = .TRUE.
            EKESIG = EKENWI
            PPRSIG = PNUCCO
            CALL SIGFER ( KPTOIP (KPROJ), EKESIG, PPRSIG, LFERMI )
            PRCOLP = ZUSFL / AUSFL * SIGMAP
            PRCOLN = ( AUSFL - ZUSFL ) / AUSFL * SIGMAN
            SIGMAA = PRCOLP + PRCOLN
            IF ( SIGMAA .GT. 0.D+00 ) THEN
               PRCOLP = PRCOLP / SIGMAA
               PRCOLN = 1.D+00 - PRCOLP
            ELSE
               PRCOLP = 0.D+00
               PRCOLN = 0.D+00
               XBIMPC = XBIMTR * TRUFAC
               YBIMPC = YBIMTR * TRUFAC
               ZBIMPC = ZBIMTR * TRUFAC
               GO TO 4500
            END IF
 5500       CONTINUE
         END IF
         XBIMPC = XBIMTR * TRUFAC
         YBIMPC = YBIMTR * TRUFAC
         ZBIMPC = ZBIMTR * TRUFAC
         IF ( COSTHE .GE. 0.D+00 ) THEN
            R0TRAJ = RSTNUC (IKPMX) * TRUFAC
            R1TRAJ = RADTOT
            IF ( R0TRAJ .GT. RADTOT ) GO TO 4500
         ELSE
            R0TRAJ = - MIN ( RSTNUC (IKPMX) * TRUFAC, RADTOT )
            R1TRAJ = RADTOT
         END IF
         CALL GRNDM(RNDM,1)
         ANMFP  = - LOG ( 1.D+00 - RNDM (1) ) / DSKRED
         IF ( BIMPCT .GT. RAD1O2 ) THEN
            SBTTSQ = 4.D+00 * ( RADTOT**2 - BIMPCT**2 ) * RHOBIM**2
            IF ( SBTTSQ .LE. ( ANMFP / SIGMAA )**2 ) GO TO 4500
         END IF
         CALL SBCOMP ( SBHAL0, SBSKI0, SBCEN0, SBCEN1, SBSKI1, SBHAL1 )
         SBTOT  = SBHAL0 + SBSKI0 + SBCEN0 + SBCEN1 + SBSKI1 + SBHAL1
         SBTOT  = RHRUSF * SBTOT
         IF ( SBTOT * SIGMAA .GT. ANMFP ) GO TO 5000
      END IF
 4500 CONTINUE
      BIMPCT = RADTOT + RADPRO + 0.0001D+00
      IF ( BIMPTR .LT. BIMPCT ) THEN
         DSTMIN = SQRT ( BIMPCT**2 - BIMPTR**2 )
         RIMPTR = BIMPCT
         XIMPTR = XBIMTR + CXIMPC * DSTMIN
         YIMPTR = YBIMTR + CYIMPC * DSTMIN
         ZIMPTR = ZBIMTR + CZIMPC * DSTMIN
      ELSE
         XIMPTR = XBIMTR
         YIMPTR = YBIMTR
         ZIMPTR = ZBIMTR
         RIMPTR = BIMPTR
      END IF
      CXIMPC = -XIMPTR / RIMPTR
      CYIMPC = -YIMPTR / RIMPTR
      CZIMPC = -ZIMPTR / RIMPTR
      RETURN
      ENTRY NWINXT ( LBCHCK )
 4800 CONTINUE
         SIGMAP = SIGMP0
         SIGMAN = SIGMN0
 4900    CONTINUE
           CALL GRNDM(RNDM,1)
           ANMFP = - LOG ( 1.D+00 - RNDM (1) ) / DSKRED
         IF ( SBRES * SIGMAA .LE. ANMFP ) GO TO 4500
 5000 CONTINUE
      SBUSED = SBUSED * RHRUSF + ANMFP / SIGMAA
      SBRES  = SBTOT  - SBUSED
      SBUSED = SBUSED / RHRUSF
      CALL GRNDM(RNDM,1)
      IF ( RNDM (1) .LT. PRCOLP ) THEN
         KNUCIM = 1
         ITFRMI = 1
      ELSE
         KNUCIM = 8
         ITFRMI = 2
      END IF
      CALL RSCOMP ( SBHAL0, SBSKI0, SBCEN0, SBCEN1, SBSKI1, SBHAL1 )
      RHOIMP = FRHONC ( ABS (RIMPCT) )
      PFRIMP = FPFRNC ( RHOIMP, ITFRMI )
      EKFIMP = FEKFNC ( PFRIMP, ITFRMI )
      RHOIMT = FRHONC ( ABS (RIMPTR) )
      PFRPRO = FPFRNC ( RHOIMT, IPWELL )
      EKFPRO = FEKFNC ( PFRPRO, IPWELL )
      VPRWLL = WLLRED * ( EKFPRO + BNDNUC )
      LELSTC = .TRUE.
      NTARGT = 1
      EKEWLL = EKENWI + VPRWLL
      EPSWLL = EKEWLL + AM (KPROJ)
      IF ( .NOT. LBCHCK ) THEN
         PPRWLL = SQRT ( EKEWLL * ( EKEWLL + 2.D+00 * AM (KPROJ) ) )
         PXPROJ = PPRWLL * CXIMPC
         PYPROJ = PPRWLL * CYIMPC
         PZPROJ = PPRWLL * CZIMPC
         PNFRMI = PFNCLV ( ITFRMI, .TRUE. )
         IF ( PNFRMI .LT. - 100.D+00 ) GO TO 4900
         CALL RACO ( PXFERM, PYFERM, PZFERM )
         PXFERM = PXFERM * PNFRMI
         PYFERM = PYFERM * PNFRMI
         PZFERM = PZFERM * PNFRMI
         EKFERM = SQRT ( PNFRMI**2 + AM (KNUCIM)**2 ) - AM (KNUCIM)
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
         EKESIG = EKEWLL
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
      RNDREJ = RNDM (1)
      IF ( RNDREJ .GE. SIGMAR ) GO TO 4800
      IF ( LBCHCK ) THEN
         ZITA   = 0.5D+00 * ( EKFIMP + EKFPRO ) / EKEWLL
         IF ( ZITA .LE. 0.5D+00 ) THEN
            PZITA = 1.D+00 - 1.4D+00 * ZITA
         ELSE
            PZITA = 1.D+00 - 1.4D+00 * ZITA + 0.4D+00 * ZITA
     &            * ( 2.D+00 - 1.D+00 / ZITA )**2.5D+00
         END IF
         RNDREJ = RNDREJ / SIGMAR
         IF ( RNDREJ .GE. PZITA ) GO TO 4800
      ELSE IF ( LFRMCK ) THEN
         EKFPRJ = EKFPRO + DEFNUC (IPWELL)
         EKFTAR = EKFIMP + DEFNUC (ITFRMI)
         IF ( EKEWLL + EKFERM .LE. EKFTAR + EKFPRJ ) GO TO 4800
         PCMSSQ = ( PXPROJ - 0.5D+00 * PXRES )**2
     &          + ( PYPROJ - 0.5D+00 * PYRES )**2
     &          + ( PZPROJ - 0.5D+00 * PZRES )**2
         RNDREJ = RNDREJ / SIGMAR
         PPROSQ = EKFPRJ * ( EKFPRJ + 2.D+00 * AM (KPROJ) )
         PTARSQ = EKFTAR * ( EKFTAR + 2.D+00 * AMNUCL (ITFRMI) )
         HLPHLP = PCMSSQ + 0.25D+00 * PTRES2
         HLPHL2 = PCMSSQ * PTRES2
         IF ( RNDREJ .GE. 0.5D+00 ) THEN
            RNDREJ = 2.D+00 * ( RNDREJ - 0.5D+00 )
            IF ( HLPHLP .LT. PTARSQ ) THEN
               COSQMX = 0.D+00
            ELSE
               COSQMX = ( HLPHLP - PTARSQ )**2 / HLPHL2
            END IF
            IF ( HLPHLP .GT. PPROSQ ) THEN
               COSQMN = 0.D+00
            ELSE
               COSQMN = ( PPROSQ - HLPHLP )**2 / HLPHL2
            END IF
         ELSE
            RNDREJ = 2.D+00 * RNDREJ
            IF ( HLPHLP .LT. PPROSQ ) THEN
               COSQMX = 0.D+00
            ELSE
               COSQMX = ( HLPHLP - PPROSQ )**2 / HLPHL2
            END IF
            IF ( HLPHLP .GT. PTARSQ ) THEN
               COSQMN = 0.D+00
            ELSE
               COSQMN = ( PTARSQ - HLPHLP )**2 / HLPHL2
            END IF
            COSQMX = ( 0.25D+00 * PTRES2 + PCMSSQ - PPROSQ )**2
     &             / ( PTRES2 * PCMSSQ )
         END IF
         RNDREJ = RNDREJ**2
         IF ( RNDREJ .GE. COSQMX .OR. RNDREJ .LT. COSQMN ) GO TO 4800
      END IF
      RETURN
*=== End of subroutine Nwisel =========================================*
      END
