*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:21:30  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.22  by  S.Giani
*-- Author :
      FUNCTION GPHSIG(Z,E)
C.
C.    ******************************************************************
C.    *                                                                *
C.    *       Function to compute the photoeffect total cross-section  *
C.    *       as a function of Z and E                                 *
C.    *           (E in GeV,  GPHSIG in barn)                          *
C.    *                                                                *
C.    *    ==>Called by : GPHOTI                                       *
C.    *       Author : L.Urban   *********                             *
C.    *                                                                *
C.    ******************************************************************
C.
#include "geant321/gconsp.inc"
      DIMENSION EK(4),EL1(4),EL2(4),C(19)
      DATA EK / 0.66644E-8 , 0.22077E-9 ,-0.32552E-11, 0.18199E-13/
      DATA EL1/-0.29179E-9 , 0.87983E-10,-0.12589E-11, 0.69602E-14/
      DATA EL2/-0.68606E-9 , 0.10078E-9 ,-0.14496E-11, 0.78809E-14/
      DATA C  /-0.88893E-6 , 0.24394E-8 , 0.28835E-6 , 0.12133E-7,
     +         -0.31104E-6 ,-0.17284E-9 , 0.14400E-7 , 0.68357E-7,
     +          0.73945E-12,-0.48149E-10, 0.55823E-9 ,-0.10089E-9,
     +         -0.10927E-5 ,-0.97897E-9 , 0.12854E-6 ,-0.45803E-5,
     +          0.16858E-10, 0.12013E-6 , 0.16924E-7/
      DATA AZ,BE/3.845,2.975/
C.
C.    ------------------------------------------------------------------
C.
      Z2=Z*Z
      X=EMASS/E
      EKZ=Z2*(EK(1)+Z*(EK(2)+Z*(EK(3)+Z*EK(4))))
      IF(E.GT.EKZ)THEN
         X2=X*X
         S=C(1)/Z+C(2)*X+C(3)+C(4)*Z+C(5)/X+C(6)*Z2+C(7)*Z/X+C(8)/X2+
     +     C(9)*Z2*Z+C(10)*Z2/X+C(11)*Z/X2+C(12)/(X2*X)
         GO TO 50
      ENDIF
C
      EL1Z=Z2*(EL1(1)+Z*(EL1(2)+Z*(EL1(3)+Z*EL1(4))))
      IF(E.GT.EL1Z)THEN
         S=C(13)/Z+C(14)*X+C(15)
         GO TO 50
      ENDIF
C
      EL2Z=Z2*(EL2(1)+Z*(EL2(2)+Z*(EL2(3)+Z*EL2(4))))
      IF(E.GT.EL2Z)THEN
         S=C(16)/Z+C(17)*X+C(18)
      ELSE
         S=C(19)
      ENDIF
C
  50  GPHSIG=S*Z**AZ*X**BE
      IF(GPHSIG.LT.0.)GPHSIG=0.
C
      END
