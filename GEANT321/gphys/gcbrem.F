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
      FUNCTION GCBREM(T)
C.    ******************************************************************
C.    *                                                                *
C.    *  It calculates the ratio of positron to electron               *
C.    *       bremsstrahlung energy loss.                              *
C.    *                                                                *
C.    *  Input:kinetic energy in GeV                                   *
C.    *  Output: ratio of positron to electron brems energy loss       *
C.    *                                                                *
C.    *  ( The ratio is calculated from a function fitted to the       *
C.    *    results of L.Kim et al.,PHYS.REV.A33,3002(1986)    )        *
C.    *                                                                *
C.    *  Author   L. Urban                                             *
C.    *                                                                *
C.    ******************************************************************
C.
#include "geant321/gcmate.inc"
#include "geant321/gconsp.inc"
      PARAMETER (AA=7522100.,A1=0.415,A3=0.0021,A5=0.00054)
C.
      CFAC=0.
      IF(Z.LE.0.) GOTO 10
      X=LOG(AA*T/Z**2)
*
      IF(X.LE.-8.) GOTO 10
*
      IF(X.GE.+9.) THEN
         CFAC=1.
      ELSE
         W=A1*X+A3*X**3+A5*X**5
         CFAC=0.5+ATAN(W)/PI
      ENDIF
*
   10 GCBREM = CFAC
*
      END
