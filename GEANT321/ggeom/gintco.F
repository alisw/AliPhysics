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
      SUBROUTINE GINTCO (X, RLEFT, RRIGHT, DZ, TAU, TAUL)
 
C        ********************************************
C        *  This subroutine finds the intersection  *
C        *  of a given ray (described by array X)   *
C        *  with a given cone (described by radii   *
C        *  RLEFT and RRIGHT and half-distance DZ). *
C        *  Output parameter is TAU, and inter-     *
C        *  section point is X = XP + TAU*XD, etc.  *
C        *                                          *
C        *  Called by GNOCON                        *
C        *  Programmed by:  Patrick Weidhaas        *
C        *          CERN,   March 1988              *
C        ********************************************
 
#include "geant321/gconsp.inc"
      DIMENSION X(6)
#if !defined(CERNLIB_SINGLE)
      DOUBLE PRECISION XP,YP,ZP,XD,YD,ZD,S,T,U,V,W
      DOUBLE PRECISION DISCR,SQDISC
#endif
C----------------------------------------------------
 
C...... Point of origin of ray:
 
      XP = X(1)
      YP = X(2)
      ZP = X(3)
 
C...... Direction cosines:
 
      XD = X(4)
      YD = X(5)
      ZD = X(6)
 
      TAU  = BIG
      TAUL = BIG
 
      S = 0.5 * (RLEFT + RRIGHT)
      T = (RLEFT - RRIGHT) / DZ
 
C......  Cone equation is:    x**2 + y**2 - Az**2 + Bz + C = 0
 
      A = 0.25 * T*T
      B = S * T
      C = -S*S
 
C......  To obtain "TAU", we must solve the quadratic equation
C......  Ut**2 + Vt + W = 0 .
 
      U = XD**2 + YD**2 - A*ZD**2
      V = 2.0 * (XP*XD + YP*YD - A*ZP*ZD) + B*ZD
      W = XP**2 + YP**2 - A*ZP**2 + B*ZP + C
 
      DISCR = V*V - 4.0*U*W
      IF (DISCR .LE. 0.0) GO TO 999
      IF(U.EQ.0.)GO TO 999
      SQDISC = SQRT (DISCR)
      TAU1 = (-V + SQDISC) / (2.0*U)
      TAU2 = (-V - SQDISC) / (2.0*U)
 
 
C......  Set TAU to the smallest positive root;
C......  otherwise let TAU = BIG .
C
C......  If both roots are positive, set TAUL to
C......  the larger one: it may be needed in the
C......  case of a PHI-segmented cone.
 
      IF (TAU1 .LT. 0.0) THEN
        IF (TAU2 .LT. 0.0) GO TO 999
        TAU = TAU2
      ELSE
        TAU = TAU1
        IF (TAU2 .GT. 0.0) THEN
          TAUL = TAU2
          IF (TAU2.LT.TAU1) THEN
            TAU = TAU2
            TAUL = TAU1
          ENDIF
        ENDIF
      ENDIF
 
  999 CONTINUE
      END
 
