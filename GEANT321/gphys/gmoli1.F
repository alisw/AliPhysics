*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:21:27  cernlib
* Geant
*
*
#include "geant321/pilot.h"
#if defined(CERNLIB_VER314)
*CMZ :  3.21/02 29/03/94  15.41.22  by  S.Giani
*-- Author :
      SUBROUTINE GMOLI1(AC,ZC,DENS,BETA2,OMC)
C.
C.    ******************************************************************
C.    *                                                                *
C.    *     Initializes material constants OMC,CHC for Moliere         *
C.    *     scattering (case of a simple material)                     *
C.    *                                                                *
C.    *    ==>Called by : GMULOF                                       *
C.    *         Author M.S. Dixit NRCC Ottawa    *********             *
C.    *                                                                *
C.    ******************************************************************
C.
C     ZETA=factor for scattering by atomic electrons
C
      DATA ZETA / 1. /
*
*     ------------------------------------------------------------------
*
      ZS=(1./AC)*ZC*(ZC+ZETA)
      ZE=-0.66666*ZS*LOG(ZC)
      ZX=ZS*LOG(1.+3.34*(ZC/(BETA2*137.036))**2)
      OMC=6702.33*ZS*DENS*EXP((ZE-ZX)/ZS)
*
      END
#endif
