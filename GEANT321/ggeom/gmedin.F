*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:20:51  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.29  by  S.Giani
*-- Author :
      SUBROUTINE GMEDIN
C.
C.    ******************************************************************
C.    *                                                                *
C.    *   Initializes common /GCVOLU/ and JGPAR structure              *
C.    *                                                                *
C.    *   Called by : GMEDIA                                           *
C.    *                                                                *
C.    *   Authors   : R.Brun, F.Bruyant                                *
C.    *                                                                *
C.    ******************************************************************
C.
#include "geant321/gcbank.inc"
#include "geant321/gcstak.inc"
#include "geant321/gcvolu.inc"
C.    ------------------------------------------------------------------
*
* *** Initialize first level in /GCVOLU/ and in JGPAR structure
*
      NLEVEL      = 1
      LVOLUM(1)   = 1
      NAMES(1)    = IQ(JVOLUM+1)
      NUMBER(1)   = 1
      LINDEX(1)   = 1
      LINMX(1)    = 1
      GONLY(1)    = 1.
      JVO         = LQ(JVOLUM-1)
      NLEVMX      = IQ(LQ(JVOLUM)-1)
      NLDEV(1)    = NLEVMX
      IF (LQ(JVO).NE.0) NLDEV(1) = 1
*
      IF (JGPAR.EQ.0) CALL MZBOOK (IXCONS, JGPAR, JGPAR, 1,'GPAR',
     +                                     NLEVMX, 0, NLEVMX, 2, 0)
      IQ(JGPAR+1) = Q(JVO+5)
      LQ(JGPAR-1) = JVO +6
*
      DO 29 I = 1,3
         GTRAN(I, 1) = 0.
         DO 19 J = 1,3
            K = (I-1)*3 +J
            GRMAT(K, 1) = 0.
   19    CONTINUE
         K = I*4 -3
         GRMAT(K, 1) = 1.
   29 CONTINUE
      GRMAT(10,1) = 0.
*                                                             END GMEDIN
      END
