*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:21:33  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/04 24/02/95  14.49.53  by  S.Giani
*-- Author :
      FUNCTION GPRSGG(Z,E)
C.
C.    ******************************************************************
C.    *                                                                *
C.    *  To calculate cross-section for e+e- production                *
C.    *  by photons (in barn/atom)                                     *
C.    *                                                                *
C.    *    ==>Called by : GPRSGA                                       *
C.    *       Author    L.Urban  *********                             *
C.    *       Modified by: G.Battistoni                                *
C.    *                                                                *
C.    *   constant cross section is assumed above 100 GeV !            *
C.    ******************************************************************
C.
#include "geant321/gconsp.inc"
      DIMENSION C(18),CC(3)
      DATA C/ 0.87842E-3,-0.19625E-2, 0.12949E-2,-0.20028E-3
     +      , 0.12575E-4,-0.28333E-6,-0.10342E-4, 0.17692E-4
     +      ,-0.82391E-5, 0.13063E-5,-0.90815E-7, 0.23586E-8
     +      ,-0.45263E-3, 0.11161E-2,-0.86749E-3, 0.21773E-3
     +      ,-0.20467E-4, 0.65372E-6/
      DATA EUPP/100./
C.
C.    ------------------------------------------------------------------
C.
      SIG=0.
*
      ESAV=E
      IF(E.GT.EUPP) E=EUPP
*
      IF(E.LE.2.*EMASS)GO TO 90
      X=LOG(E/EMASS)
      DO 20 I=1,3
         CC(I)=0.
         D=1.
         DO 10 J=1,6
            JJ=6*I+J-6
            CC(I)=CC(I)+D*C(JJ)
            D=D*X
  10     CONTINUE
  20  CONTINUE
C
      SIG=Z*(Z+1.)*(CC(1)+CC(2)*Z+CC(3)/Z)
C
  90  IF(SIG.LT.0.)SIG=0.
      GPRSGG=SIG
*
      E=ESAV
*
      END
