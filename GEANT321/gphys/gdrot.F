*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:21:24  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.21  by  S.Giani
*-- Author :
      SUBROUTINE GDROT(P,COSTH,SINTH,COSPH,SINPH)
C.
C.    ******************************************************************
C.    *                                                                *
C.    *  Rotates vector from one reference system into another.        *
C.    *  THETA and PHI are anti-clockwise Eulerian angles between the  *
C.    *  two systems.                                                  *
C.    *                                                                *
C.    *    ==>Called by : GBREME,GCOMP,GDECAY,GDECA3,GDRAY,GPAIRG      *
C.    *       Author    M.Hansroul, G.Patrick  *********               *
C.    *                                                                *
C.    *                                                                *
C.    ******************************************************************
C.
      DIMENSION P(3)
#if !defined(CERNLIB_SINGLE)
      DOUBLE PRECISION P1,P2,P3
#endif
C.
C.    ------------------------------------------------------------------
C.
      P1=P(1)
      P2=P(2)
      P3=P(3)
      P(1)=P1*COSTH*COSPH - P2*SINPH + P3*SINTH*COSPH
      P(2)=P1*COSTH*SINPH + P2*COSPH + P3*SINTH*SINPH
      P(3)=-P1*SINTH                 + P3*COSTH
C
      END
