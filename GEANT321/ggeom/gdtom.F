*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:20:48  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.28  by  S.Giani
*-- Author :
      SUBROUTINE GDTOM(XD,XM,IFLAG)
C.
C.    ******************************************************************
C.    *                                                                *
C.    *       Computes coordinates XM (Master Reference System         *
C.    *       knowing the coordinates XD (Detector Ref System)         *
C.    *       The local reference system can be initialized by         *
C.    *         - the tracking routines and GDTOM used in GUSTEP       *
C.    *         - a call to GSCMED(NLEVEL,NAMES,NUMBER)                *
C.    *             (inverse routine is GMTOD)                         *
C.    *                                                                *
C.    *        If IFLAG=1  convert coordinates                         *
C.    *           IFLAG=2  convert direction cosinus                   *
C.    *                                                                *
C.    *    ==>Called by : <USER>                                       *
C.    *       Author    R.Brun ,M Hansroul  *********                  *
C.    *                                                                *
C.    ******************************************************************
C.
#include "geant321/gcvolu.inc"
      DIMENSION XD(3),XM(3)
      EQUIVALENCE (N,NLEVEL)
C.
C.    ------------------------------------------------------------------
C.
      IF(IFLAG.EQ.1)THEN
         IF(GRMAT(10,N).NE.0.)THEN
            XM(1) = GTRAN(1,N) + GRMAT(1,N)*XD(1) + GRMAT(4,N)*XD(2) +
     +      GRMAT(7 ,N)*XD(3)
            XM(2) = GTRAN(2,N) + GRMAT(2,N)*XD(1) + GRMAT(5,N)*XD(2) +
     +      GRMAT(8 ,N)*XD(3)
            XM(3) = GTRAN(3,N) + GRMAT(3,N)*XD(1) + GRMAT(6,N)*XD(2) +
     +      GRMAT(9 ,N)*XD(3)
*
         ELSE
            XM(1) = GTRAN(1,N) + XD(1)
            XM(2) = GTRAN(2,N) + XD(2)
            XM(3) = GTRAN(3,N) + XD(3)
*
         ENDIF
*
      ELSE
         IF(GRMAT(10,N).NE.0.)THEN
            XM(1)=GRMAT(1,N)*XD(1)+GRMAT(4,N)*XD(2)+GRMAT(7,N)*XD(3)
            XM(2)=GRMAT(2,N)*XD(1)+GRMAT(5,N)*XD(2)+GRMAT(8,N)*XD(3)
            XM(3)=GRMAT(3,N)*XD(1)+GRMAT(6,N)*XD(2)+GRMAT(9,N)*XD(3)
*
         ELSE
            XM(1) = XD(1)
            XM(2) = XD(2)
            XM(3) = XD(3)
*
         ENDIF
*
      ENDIF
      END
