*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:22:02  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.46  by  S.Giani
*-- Author :
*$ CREATE PHDWLL.FOR
*COPY PHDWLL
*
*=== phdwll ===========================================================*
*
      SUBROUTINE PHDWLL ( UBIMPT, VBIMPT, WBIMPT )
 
#include "geant321/dblprc.inc"
#include "geant321/dimpar.inc"
#include "geant321/iounit.inc"
*
*----------------------------------------------------------------------*
*----------------------------------------------------------------------*
*
#include "geant321/nucgeo.inc"
#include "geant321/parnuc.inc"
#include "geant321/part.inc"
*
 1000 CONTINUE
      PDIFF = PNUCCO - PPRWLL
      IF ( PDIFF .LT. - ANGLGB ) THEN
         IF ( RIMPTR .LE. RADIU0 ) THEN
            RADHLP = 0.5D+00 * ( RADTOT + RADPRO + MAX ( ABS (RIMPTR),
     &               RADIU0 ) )
            CZHLP  = SQRT ( ( RADHLP + BIMPTR ) * ( RADHLP - BIMPTR ) )
     &             / RADHLP
            HLPHLP = RIMPTR / ( RIMPCT * RADHLP )
            CXHLP  = CZHLP * CXIMPC - XBIMPC * HLPHLP
            CYHLP  = CZHLP * CYIMPC - YBIMPC * HLPHLP
            CZHLP  = CZHLP * CZIMPC - ZBIMPC * HLPHLP
            PXPROJ = PNUCCO * CXIMPC
            PYPROJ = PNUCCO * CYIMPC
            PZPROJ = PNUCCO * CZIMPC
            PDTCMP = PXPROJ * CXHLP + PYPROJ * CYHLP + PZPROJ * CZHLP
            DELTAE = PDTCMP**2 - PNUCCO**2 + PPRWLL**2
            DELTAP = - PDTCMP + SQRT ( DELTAE )
            PXPROJ = PXPROJ + DELTAP * CXHLP
            PYPROJ = PYPROJ + DELTAP * CYHLP
            PZPROJ = PZPROJ + DELTAP * CZHLP
         ELSE
            EKEBIM = MAX ( EKECON + VPRBIM, EKEWLL )
            PBIMSQ = EKEBIM * ( EKEBIM + 2.D+00 * AM (KPRIN) )
            RADHLP = 0.5D+00 * ( RADTOT + RADPRO + MAX ( BIMPTR,
     &               RADIU0 ) )
            CZHLP  = SQRT ( ( RADHLP + BIMPTR ) * ( RADHLP - BIMPTR ) )
     &             / RADHLP
            HLPHLP = RIMPTR / ( RIMPCT * RADHLP )
            CXHLP  = CZHLP * CXIMPC - XBIMPC * HLPHLP
            CYHLP  = CZHLP * CYIMPC - YBIMPC * HLPHLP
            CZHLP  = CZHLP * CZIMPC - ZBIMPC * HLPHLP
            PXPROJ = PNUCCO * CXIMPC
            PYPROJ = PNUCCO * CYIMPC
            PZPROJ = PNUCCO * CZIMPC
            PDTCMP = PXPROJ * CXHLP + PYPROJ * CYHLP + PZPROJ * CZHLP
            DELTAE = PDTCMP**2 - PNUCCO**2 + PBIMSQ
            DELTAP = - PDTCMP + SQRT ( DELTAE )
            PXPROJ = PXPROJ + DELTAP * CXHLP
            PYPROJ = PYPROJ + DELTAP * CYHLP
            PZPROJ = PZPROJ + DELTAP * CZHLP
            PPBIM  = SQRT ( PBIMSQ )
            COSTHE = ( PXPROJ * CXIMPC + PYPROJ * CYIMPC
     &             + PZPROJ * CZIMPC ) / PPBIM
            THETA  = ACOS (COSTHE) * ( 1.D+00 + ( PNUCCO - PPBIM )
     &             / PDIFF )
            SINTHE = SIN (THETA)
            COSTHE = COS (THETA)
            PXPROJ = PPRWLL * ( COSTHE * CXIMPC - SINTHE * UBIMPT )
            PYPROJ = PPRWLL * ( COSTHE * CYIMPC - SINTHE * VBIMPT )
            PZPROJ = PPRWLL * ( COSTHE * CZIMPC - SINTHE * WBIMPT )
         END IF
      ELSE IF ( PDIFF .GT. ANGLGB ) THEN
         STOP 'PHDWLL'
      ELSE
         PXPROJ = PPRWLL * CXIMPC
         PYPROJ = PPRWLL * CYIMPC
         PZPROJ = PPRWLL * CZIMPC
      END IF
      RETURN
*=== End of subroutine phdwll =========================================*
      END
