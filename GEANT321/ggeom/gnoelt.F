*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:20:52  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.29  by  S.Giani
*-- Author :
      SUBROUTINE GNOELT(X,PAR,IACT,SNEXT,SNXT,SAFE)
C
C     ****************************************************************
C     *                                                              *
C     *     Compute distance up to intersection with 'ELTU' volume,  *
C     *     from outside point X(1-3) along direction X(4-6).        *
C     *                                                              *
C     *     PAR    (input)  : volume parameters                      *
C     *     IACT   (input)  : action flag                            *
C     *       = 0   Compute SAFE only                                *
C     *       = 1   Compute SAFE, and SNXT only if SNEXT.gt.new SAFE *
C     *       = 2   Compute both SAFE and SNXT                       *
C     *       = 3   Compute SNXT only                                *
C     *     SNEXT  (input)  : see IACT = 1                           *
C     *     SNXT   (output) : distance to volume boundary            *
C     *     SAFE   (output) : shortest distance to any boundary      *
C     *                                                              *
C     *  ==>Called by : GNEXT,GTNEXT                                 *
C     *       Author  A.Solano                                       *
C     *                                                              *
C     ****************************************************************
C
#include "geant321/gconsp.inc"
C
      DIMENSION X(6),PAR(3),TAU(2)
C
#if !defined(CERNLIB_SINGLE)
      DOUBLE PRECISION SAFZ,A2,B2,X0,Y0,PHI1,PHI2,PHI3,X3,Y3
      DOUBLE PRECISION DXY2,U,V,W,DISCR,SQDISC,TAU,TAUZ,ZI,XZ,YZ
#endif
      SNXT = BIG
      SAFZ = ABS(X(3))-PAR(3)
      A2 = PAR(1)*PAR(1)
      B2 = PAR(2)*PAR(2)
C
      IF(IACT.EQ.3)GOTO 40
C
C      -----------------------------------
C      |  Compute safety-distance 'SAFE' |
C      -----------------------------------
C
C ....  First check Z
      X0 = ABS(X(1))
      Y0 = ABS(X(2))
C
      SAFE=0.
      IF(X0*X0/A2+Y0*Y0/B2.LT.1.) GO TO 30
      PHI1=0.
      PHI2=PIBY2
      DO 10    I=1,10
         PHI3=(PHI1+PHI2)*0.5
         X3=PAR(1)*COS(PHI3)
         Y3=PAR(2)*SIN(PHI3)
         D=Y3*A2*(X0-X3)-X3*B2*(Y0-Y3)
*
*        D is the (signed) distance from the point (X0,Y0)
*        to the normal to the ellipse at the point (X3,Y3).
*
         IF (D.LT.0.) THEN
            PHI1=PHI3
         ELSE
            PHI2=PHI3
         END IF
   10 CONTINUE
   20 SAFE=SQRT((X0-X3)**2+(Y0-Y3)**2)-.01
   30 IF(SAFZ.GT.0.)THEN
*
* ....   Combine the radial distance whit the Z-distance
*
         SAFE = SQRT(SAFE**2+SAFZ**2)
      ENDIF
      IF(IACT.EQ.0)GOTO 999
      IF(IACT.EQ.1.AND.SNEXT.LT.SAFE)GOTO 999
   40 CONTINUE
C
C      ---------------------------------------
C      |  Compute the vector-distance 'SNXT' |
C      ---------------------------------------
C
      IF(SAFZ.GT.0.0.AND.X(3)*X(6).GE.0.0)GOTO 999
C
      DXY2 = (1-X(6))*(1+X(6))
      IF(DXY2.LE.0.)GOTO 60
C
C ....  Find the intersection of the given ray
C       (described by array X) whit the cylider.
C       Ray equation : X'(1-3) = X(1-3) + TAU*X(4-6)
C       Cylinder equation : x**2/a**2 + y**2/b**2 = 1
C       To obtain TAU,solve the quadratic equation
C       Ut**2 + 2Vt + W = 0
C
      U = X(4)*X(4)*B2+X(5)*X(5)*A2
      V = X(1)*X(4)*B2+X(2)*X(5)*A2
      W = X(1)*X(1)*B2+X(2)*X(2)*A2-A2*B2
      DISCR = V*V-U*W
      IF(DISCR.LT.0.)GOTO 999
      IF(U.EQ.0.)GOTO 999
      SQDISC = SQRT(DISCR)
      TAU(1) = (-V+SQDISC)/U
      TAU(2) = (-V-SQDISC)/U
C
      DO 50 I=1,2
         IF(TAU(I).GE.0.)THEN
            ZI = X(3)+TAU(I)*X(6)
            IF((ABS(ZI)-PAR(3)).LT.1.E-6)THEN
C
C ....  Set SNXT to the smallest positive TAU,only if
C       the intersection point is inside the Z limits
C
               SNXT = MIN(SNXT,REAL(TAU(I)))
            ENDIF
         ENDIF
   50 CONTINUE
   60 CONTINUE
C
      IF(SAFZ.GT.0.)THEN
C
C ....  Check intersection whit Z planes
C
         IF(X(3).GT.0.) ZI = PAR(3)
         IF(X(3).LT.0.) ZI = -PAR(3)
C
         TAUZ = (ZI-X(3))/X(6)
         XZ = X(1)+X(4)*TAUZ
         YZ = X(2)+X(5)*TAUZ
         IF((XZ*XZ/A2+YZ*YZ/B2).LE.1.) SNXT = TAUZ
      ENDIF
C
  999 END
