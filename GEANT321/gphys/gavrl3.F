*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:21:21  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.21  by  S.Giani
*-- Author :
      REAL FUNCTION GAVRL3(GAMA,BETA,XM)
C.
C.    ******************************************************************
C.    *                                                                *
C.    *  L3 shell x-section a la Gavrila                               *
C.    *    Gama  - L3 shell photoelectron gama                         *
C.    *    Beta  - L3 shell photoelectron beta                         *
C.    *    XM    - ratio of the electron mass to the photon energy     *
C.    *                                                                *
C.    *    ==>CALLED BY : GPHOT                                        *
C.    *       AUTHOR    : J. Chwastowski                               *
C.    *                                                                *
C.    ******************************************************************
C.
#if !defined(CERNLIB_SINGLE)
      DOUBLE PRECISION G,G2,G3,B,P1,P2,P3
#endif
      G = GAMA
      G2 = G*G
      G3 = G2*G
      B  = BETA
      P1 = (4*G3-6*G2+5*G+3)
      P2 = -(G2-3*G+4)
      P3 = 1
      IF(BETA.GT.0.0001) P3 = LOG(G*(1+B))/(B*G)
      GAVRL3 = (P1+P2*P3)*B*G*XM**5
      END
