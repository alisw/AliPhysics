*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:21:23  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.21  by  S.Giani
*-- Author :
      SUBROUTINE GCOMPI
C.
C.    ******************************************************************
C.    *                                                                *
C.    *  Calculates cross-section in current material                  *
C.    *  for COMPTON scattering.                                       *
C.    *                                                                *
C.    *  The total cross-section of L.Urban is used to estimate the    *
C.    *  photon mean free path.                                        *
C.    *    It gives a good description of the data from 10KeV to       *
C.    *    100/Z GeV                                                   *
C.    *                                                                *
C.    *    ==>Called by : GPHYSI                                       *
C.    *       Authors    R.Brun, L.Urban  *********                    *
C.    *                                                                *
C.    ******************************************************************
C.
#include "geant321/gcbank.inc"
#include "geant321/gctrak.inc"
#include "geant321/gconsp.inc"
#include "geant321/gcmate.inc"
#include "geant321/gcjloc.inc"
#include "geant321/gcmulo.inc"
      DIMENSION PCOMP(12),P(4)
      DATA AC,BC,CC/20.,230.,440./
      DATA PCOMP/ 0.27965,-0.18300,0.67527E+1,-0.19798E+2,
     +            0.19756E-4,-0.10205E-1,-0.73913E-1,0.27079E-1,
     +           -0.39178E-6,0.68241E-4,0.60480E-4,0.30274E-3/
      DATA PI2R2I/2.004275/
C.
C.    ------------------------------------------------------------------
C.
C             Energy cut-off above which we ignore COMPTON scattering.
C
      SIG=0.
      IF(Z.LT.1.)GO TO 90
      EGAM   = ELOW(IEKBIN)
      X=EGAM/EMASS
      IF(EGAM.GT.(100./Z))GO TO 90
      X2=X*X
      JPROB=LQ(JMA-4)
C
      DO 20 I=1,4
         C=1.
         P(I)=0.
         DO 10 J=1,3
            JJ=4*J+I-4
            P(I)=P(I)+C*PCOMP(JJ)
            C=C*Z
  10     CONTINUE
  20  CONTINUE
C
      SIG=P(1)*LOG(1.+2.*X)/X
      SIG=SIG+(P(2)+P(3)*X+P(4)*X2)/(1.+AC*X+BC*X2+CC*X2*X)
      SIG=PI2R2I*Q(JPROB+17)*SIG
C
  90  IF(SIG.GT.0.)THEN
         Q(JCOMP+IEKBIN)=1./SIG
      ELSE
         Q(JCOMP+IEKBIN)=BIG
      ENDIF
C
      END
