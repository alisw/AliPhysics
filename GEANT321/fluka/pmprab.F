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
*$ CREATE PMPRAB.FOR
*COPY PMPRAB
*
*=== pmprab ===========================================================*
*
      SUBROUTINE PMPRAB ( KPROJ, EKIN, PPROJ, TXX, TYY, TZZ, WEE )
 
#include "geant321/dblprc.inc"
#include "geant321/dimpar.inc"
#include "geant321/iounit.inc"
*
*----------------------------------------------------------------------*
*                                                                      *
*     Created on 22 september 1991 by    Alfredo Ferrari & Paola Sala  *
*                                                   Infn - Milan       *
*                                                                      *
*     Last change on 22-sep-91     by    Alfredo Ferrari               *
*                                                                      *
*                                                                      *
*----------------------------------------------------------------------*
*
#include "geant321/balanc.inc"
#include "geant321/finuc.inc"
#include "geant321/nucdat.inc"
#include "geant321/nucgeo.inc"
#include "geant321/parevt.inc"
#include "geant321/parnuc.inc"
#include "geant321/part.inc"
#include "geant321/resnuc.inc"
*
      REAL RNDM(1)
*
      IF ( KPROJ .NE. 14 .OR. EKIN .GT. 2.D+00 * GAMMIN .OR. IBTAR .NE.
     &     1 .OR. ICHTAR .NE. 1 ) THEN
         WRITE (LUNOUT,*)' **** Pmprab: kproj,ekin,ibtar,ichtar',
     &                     KPROJ,EKIN,IBTAR,ICHTAR
         WRITE (LUNERR,*)' **** Pmprab: kproj,ekin,ibtar,ichtar',
     &                     KPROJ,EKIN,IBTAR,ICHTAR
      END IF
      PXRES  = PXTTOT
      PYRES  = PYTTOT
      PZRES  = PZTTOT
      PTRES  = PTTOT
      CALL GRNDM(RNDM,1)
      RNDPAN = RNDM (1)
      IF ( RNDPAN .GE. 1.D+00 / PNFRAT ) THEN
         ERES   = EKIN + AM (KPROJ) + EKFERM + AM (1)
         UMO2   = ( ERES - PTRES ) * ( ERES + PTRES )
         UMO    = SQRT (UMO2)
         GAMCM = ERES  / UMO
         ETAX  = PXRES / UMO
         ETAY  = PYRES / UMO
         ETAZ  = PZRES / UMO
         ECMSNU = 0.5D+00 * ( UMO2 + AMNUSQ (2) ) / UMO
         PCMS   = UMO - ECMSNU
         CALL RACO ( PCMSX, PCMSY, PCMSZ )
         PCMSX = PCMS * PCMSX
         PCMSY = PCMS * PCMSY
         PCMSZ = PCMS * PCMSZ
         ETAPCM = ETAX * PCMSX + ETAY * PCMSY + ETAZ * PCMSZ
         PHELP  = PCMS + ETAPCM / ( GAMCM + 1.D+00 )
         ENPHOT = GAMCM * PCMS + ETAPCM
         PXHELP = PCMSX + ETAX * PHELP
         PYHELP = PCMSY + ETAY * PHELP
         PZHELP = PCMSZ + ETAZ * PHELP
         PXRES = PXRES - PXHELP
         PYRES = PYRES - PYHELP
         PZRES = PZRES - PZHELP
         ERES  = ERES  - ENPHOT
         NP = NP + 1
         TKI   (NP) = ENPHOT
         KPART (NP) = 7
         PLR   (NP) = ENPHOT
         CXR   (NP) = PXHELP / PLR (NP)
         CYR   (NP) = PYHELP / PLR (NP)
         CZR   (NP) = PZHELP / PLR (NP)
         WEI   (NP) = WEE
         IOTHER = IOTHER + 1
         PXNUCR = PXNUCR + PXHELP
         PYNUCR = PYNUCR + PYHELP
         PZNUCR = PZNUCR + PZHELP
         ENUCR  = ENUCR  + TKI (NP)
         IBNUCR = IBNUCR + IBAR (KPART(NP))
         ICNUCR = ICNUCR + ICH  (KPART(NP))
         ETAPCM = - ETAPCM
         PHELP  = ECMSNU + ETAPCM / ( GAMCM + 1.D+00 )
         ENNEU  = GAMCM * ECMSNU + ETAPCM
         PXHELP = -PCMSX + ETAX * PHELP
         PYHELP = -PCMSY + ETAY * PHELP
         PZHELP = -PCMSZ + ETAZ * PHELP
         NP = NP + 1
         TKI   (NP) = ENNEU - AM (8)
         KPART (NP) = 8
         PLR (NP) = SQRT ( PXHELP**2 + PYHELP**2 + PZHELP**2 )
         CXR   (NP) = PXHELP / PLR (NP)
         CYR   (NP) = PYHELP / PLR (NP)
         CZR   (NP) = PZHELP / PLR (NP)
         WEI   (NP) = WEE
         ERES  = ERES  - ENNEU
         PXRES = PXRES - PXHELP
         PYRES = PYRES - PYHELP
         PZRES = PZRES - PZHELP
         IBRES = 0
         ICRES = 0
         PTRES = 0.D+00
         ANOW  = 0.D+00
         ZNOW  = 0.D+00
      ELSE
         ERES   = EKIN + AM (KPROJ) + EKFERM + AM (1)
         UMO2   = ( ERES - PTRES ) * ( ERES + PTRES )
         UMO    = SQRT (UMO2)
         GAMCM = ERES  / UMO
         ETAX  = PXRES / UMO
         ETAY  = PYRES / UMO
         ETAZ  = PZRES / UMO
         ECMSNU = 0.5D+00 * ( UMO2 + AM (8)* AM (8) - AM (23) * AM (23)
     &          ) / UMO
         ECMSP0 = UMO - ECMSNU
         PCMS = SQRT ( ( ECMSP0 - AM (23) ) * ( ECMSP0 + AM (23) ) )
         CALL RACO ( PCMSX, PCMSY, PCMSZ )
         PCMSX = PCMS * PCMSX
         PCMSY = PCMS * PCMSY
         PCMSZ = PCMS * PCMSZ
         ETAPCM = ETAX * PCMSX + ETAY * PCMSY + ETAZ * PCMSZ
         PHELP  = ECMSP0 + ETAPCM / ( GAMCM + 1.D+00 )
         ENPIO0 = GAMCM * ECMSP0 + ETAPCM
         PXHELP = PCMSX + ETAX * PHELP
         PYHELP = PCMSY + ETAY * PHELP
         PZHELP = PCMSZ + ETAZ * PHELP
         PXRES = PXRES - PXHELP
         PYRES = PYRES - PYHELP
         PZRES = PZRES - PZHELP
         ERES  = ERES  - ENPIO0
         NP = NP + 1
         TKI   (NP) = ENPIO0 - AM (23)
         KPART (NP) = 23
         PLR (NP) = SQRT ( PXHELP**2 + PYHELP**2 + PZHELP**2 )
         CXR   (NP) = PXHELP / PLR (NP)
         CYR   (NP) = PYHELP / PLR (NP)
         CZR   (NP) = PZHELP / PLR (NP)
         WEI   (NP) = WEE
         IOTHER = IOTHER + 1
         PXNUCR = PXNUCR + PXHELP
         PYNUCR = PYNUCR + PYHELP
         PZNUCR = PZNUCR + PZHELP
         ENUCR  = ENUCR  + TKI (NP)
         IBNUCR = IBNUCR + IBAR (KPART(NP))
         ICNUCR = ICNUCR + ICH  (KPART(NP))
         ETAPCM = - ETAPCM
         PHELP  = ECMSNU + ETAPCM / ( GAMCM + 1.D+00 )
         ENNEU  = GAMCM * ECMSNU + ETAPCM
         PXHELP = -PCMSX + ETAX * PHELP
         PYHELP = -PCMSY + ETAY * PHELP
         PZHELP = -PCMSZ + ETAZ * PHELP
         NP = NP + 1
         TKI   (NP) = ENNEU - AM (8)
         KPART (NP) = 8
         PLR (NP) = SQRT ( PXHELP**2 + PYHELP**2 + PZHELP**2 )
         CXR   (NP) = PXHELP / PLR (NP)
         CYR   (NP) = PYHELP / PLR (NP)
         CZR   (NP) = PZHELP / PLR (NP)
         WEI   (NP) = WEE
         ERES  = ERES  - ENNEU
         PXRES = PXRES - PXHELP
         PYRES = PYRES - PYHELP
         PZRES = PZRES - PZHELP
         IBRES = 0
         ICRES = 0
         PTRES = 0.D+00
         ANOW  = 0.D+00
         ZNOW  = 0.D+00
      END IF
      RETURN
*=== End of subroutine PMPRAB =========================================*
      END
