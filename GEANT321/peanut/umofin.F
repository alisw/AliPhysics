*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:22:04  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.46  by  S.Giani
*-- Author :
*$ CREATE UMOFIN.FOR
*COPY UMOFIN
*
*=== umofin ===========================================================*
*
      SUBROUTINE UMOFIN ( IKPMX, BBRES, ZZRES, LTRPPD )
 
#include "geant321/dblprc.inc"
#include "geant321/dimpar.inc"
#include "geant321/iounit.inc"
*
*----------------------------------------------------------------------*
*----------------------------------------------------------------------*
*
#include "geant321/balanc.inc"
#include "geant321/nucdat.inc"
#include "geant321/nucgeo.inc"
#include "geant321/parevt.inc"
#include "geant321/parnuc.inc"
#include "geant321/part.inc"
#include "geant321/resnuc.inc"
*
      LOGICAL LTRPPD
*
      PXRES = PXTTOT - PXNUCR - PXINTR - PXNUCL (IKPMX)
      PYRES = PYTTOT - PYNUCR - PYINTR - PYNUCL (IKPMX)
      PZRES = PZTTOT - PZNUCR - PZINTR - PZNUCL (IKPMX)
      ERES  = ETTOT  - EINTR  - ENUCR  - ENNUC  (IKPMX)
      PTRES2 = PXRES**2 + PYRES**2 + PZRES**2
      UMO2   = ERES**2 - PTRES2
      DELTU2 = ( AMNRES + EEXMIN )**2 - UMO2
      IF ( DELTU2 .LT. 0.D+00 ) RETURN
*
      PNUCC0 = PNUCCO
      ENNUC0 = ENNUC (IKPMX)
*
      PXRES = PXTTOT - PXNUCR - PXINTR
      PYRES = PYTTOT - PYNUCR - PYINTR
      PZRES = PZTTOT - PZNUCR - PZINTR
      ERES  = ETTOT  - EINTR  - ENUCR
      PTRES2 = PXRES**2 + PYRES**2 + PZRES**2
      UMO2   = ERES**2  - PTRES2
      AMSQ   = AM (KPRIN)**2
      ENEMIN = 0.5D+00 * ( UMO2 - ( AMNRES + EEXMIN )**2 - AMSQ )
     &       / ( AMNRES + EEXMIN )
      ENEDLT = 0.5D+00 * ( UMO2 - ( AMNRES + EEXDEL )**2 - AMSQ )
     &       / ( AMNRES + EEXDEL )
      ENEANY = 0.5D+00 * ( UMO2 - ( AMNRES + EEXANY )**2 - AMSQ )
     &       / ( AMNRES + EEXANY )
      IF ( ENEANY .LE. AM (KPRIN) ) THEN
         IBRES = IBRES + IBAR (KPRIN)
         ICRES = ICRES + ICH  (KPRIN)
         BBRES = IBRES
         ZZRES = ICRES
         AMMRES = BBRES * AMUAMU + 0.001D+00
     &          * FKENER ( BBRES, ZZRES )
         AMNRES = AMMRES - ZZRES * AMELEC + ELBNDE ( ICRES )
         LTRPPD = .TRUE.
         EEXMIN = 0.D+00
         RETURN
      END IF
      IF ( ENNUC (IKPMX) .GT. 0.5D+00 * ( ENEANY + ENEDLT ) ) THEN
         ENNUC (IKPMX) = ENEANY
         PNUCL (IKPMX) = SQRT ( ( ENNUC (IKPMX) - AM (KPRIN) )
     &                 * ( ENNUC (IKPMX) + AM (KPRIN) ) )
         PHELP = PNUCL (IKPMX) / PNUCC0
         PXNUCL (IKPMX) = PXNUCL (IKPMX) * PHELP
         PYNUCL (IKPMX) = PYNUCL (IKPMX) * PHELP
         PZNUCL (IKPMX) = PZNUCL (IKPMX) * PHELP
      ELSE IF ( ENNUC (IKPMX) .GT. 0.5D+00 * ( ENEDLT + ENEMIN ) ) THEN
         ENNUC (IKPMX) = ENEDLT
         PNUCL (IKPMX) = SQRT ( ( ENNUC (IKPMX) - AM (KPRIN) )
     &                 * ( ENNUC (IKPMX) + AM (KPRIN) ) )
         PHELP = PNUCL (IKPMX) / PNUCC0
         PXNUCL (IKPMX) = PXNUCL (IKPMX) * PHELP
         PYNUCL (IKPMX) = PYNUCL (IKPMX) * PHELP
         PZNUCL (IKPMX) = PZNUCL (IKPMX) * PHELP
      ELSE IF ( ENNUC (IKPMX) .GE. ENEMIN ) THEN
         ENNUC (IKPMX) = ENEMIN
         PNUCL (IKPMX) = SQRT ( ( ENNUC (IKPMX) - AM (KPRIN) )
     &                 * ( ENNUC (IKPMX) + AM (KPRIN) ) )
         PHELP = PNUCL (IKPMX) / PNUCC0
         PXNUCL (IKPMX) = PXNUCL (IKPMX) * PHELP
         PYNUCL (IKPMX) = PYNUCL (IKPMX) * PHELP
         PZNUCL (IKPMX) = PZNUCL (IKPMX) * PHELP
      ELSE
      END IF
      UMO = SQRT (UMO2)
      AMNHLP = SQRT ( UMO2 + PNUCL (IKPMX)**2 ) - ENNUC (IKPMX)
      GAMLAB = ERES  / UMO
      ETXLAB = PXRES / UMO
      ETYLAB = PYRES / UMO
      ETZLAB = PZRES / UMO
      GAMRES = ( AMNHLP + ENNUC (IKPMX) ) / UMO
      ETXRES = PXNUCL (IKPMX) / UMO
      ETYRES = PYNUCL (IKPMX) / UMO
      ETZRES = PZNUCL (IKPMX) / UMO
      ETAPCM = ETXRES * PXNUCL (IKPMX) + ETYRES * PYNUCL (IKPMX)
     &       + ETZRES * PZNUCL (IKPMX)
      ECMSPR = GAMRES * ENNUC (IKPMX) - ETAPCM
      PHELP  = ENNUC (IKPMX) - ETAPCM / ( GAMRES + 1.D+00 )
      PCMSX  = PXNUCL (IKPMX) - ETXRES * PHELP
      PCMSY  = PYNUCL (IKPMX) - ETYRES * PHELP
      PCMSZ  = PZNUCL (IKPMX) - ETZRES * PHELP
      ETAPCM = ETXLAB * PCMSX + ETYLAB * PCMSY + ETZLAB * PCMSZ
      ENNUC  (IKPMX) = GAMLAB * ECMSPR + ETAPCM
      PHELP  = ECMSPR + ETAPCM / ( GAMLAB + 1.D+00 )
      PXNUCL (IKPMX) = PCMSX + ETXLAB * PHELP
      PYNUCL (IKPMX) = PCMSY + ETYLAB * PHELP
      PZNUCL (IKPMX) = PCMSZ + ETZLAB * PHELP
      PNUCL  (IKPMX) = SQRT ( ( ENNUC (IKPMX) - AM (KPRIN) )
     &               * ( ENNUC (IKPMX) + AM (KPRIN) ) )
 
      RETURN
*=== End of subroutine UMOFIN =========================================*
      END
