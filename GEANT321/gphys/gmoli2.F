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
      SUBROUTINE GMOLI2(AC,ZC,WMAT,NLM,DENS,BETA2,OMC)
C.
C.    ******************************************************************
C.    *                                                                *
C.    *     Initializes material constants OMC,CHC for Moliere         *
C.    *     scattering.                                                *
C.    *                                                                *
C.    * WMAT=proportion by weight of NLM  individual constituents of   *
C.    *      charge ZC and atomic weight AC for the mixture/compound   *
C.    * DENS=density in gm/cm**3                                       *
C.    *                                                                *
C.    *    ==>Called by : GMULOF                                       *
C.    *         Author M.S. Dixit NRCC Ottawa    *********             *
C.    *                                                                *
C.    ******************************************************************
C.
C     ZETA=factor for scattering by atomic electrons
C
      DIMENSION WMAT(1),AC(1),ZC(1)
      DATA ZETA / 1. /
*
*     ------------------------------------------------------------------
*
      ZS=0.
      ZE=0.
      ZX=0.
      DO 15 I=1,NLM
         TMP=(WMAT(I)/AC(I))*ZC(I)*(ZC(I)+ZETA)
         ZS=ZS+TMP
         ZE=ZE-0.66666*TMP*LOG(ZC(I))
         ZX=ZX+TMP*LOG(1.+3.34*(ZC(I)/(BETA2*137.036))**2)
15    CONTINUE
      OMC=6702.33*ZS*DENS*EXP((ZE-ZX)/ZS)
*
      END
 
#endif
