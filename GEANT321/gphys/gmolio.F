*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:21:27  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.22  by  S.Giani
*-- Author :
      SUBROUTINE GMOLIO(AC,ZC,WMAT,NLM,DENS,BETA2,CHARG2,OMC)
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
      DIMENSION WMAT(*),AC(*),ZC(*)
      PARAMETER (ZETA=1.)
*
* *** This is equal to 3.34/137.036**2
*
      PARAMETER (CONST1 = 1.778595E-4)
*
*     ------------------------------------------------------------------
*
      ZS     = 0.
      ZE     = 0.
      ZX     = 0.
      CONST  = CONST1*CHARG2/BETA2
      DO 10 I=1,NLM
         TMP   = (WMAT(I)/AC(I))*ZC(I)*(ZC(I)+ZETA)
         ZS    = ZS+TMP
         ZE    = ZE-0.66666*TMP*LOG(ZC(I))
         ZX    = ZX+TMP*LOG(1.+CONST*ZC(I)**2)
  10  CONTINUE
      OMC    = 6702.33*ZS*DENS*EXP((ZE-ZX)/ZS) * CHARG2
*
      END
 
