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
      SUBROUTINE GNELTU(X,PAR,IACT,SNEXT,SNXT,SAFE)
C
C     ****************************************************************
C     *                                                              *
C     *     Compute distance up to intersection with 'ELTU' volume,  *
C     *     from inside point X(1-3) along direction X(4-6).         *
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
      DIMENSION X(6),PAR(3)
#if !defined(CERNLIB_SINGLE)
      DOUBLE PRECISION SAFZ1,SAFZ2,SAFR,A2,B2,X0,Y0,A,B,D1,D2
      DOUBLE PRECISION X1,X2,X3,Y1,Y2,Y3
      DOUBLE PRECISION SZ,XZ,YZ,U,V,W,DISCR,SQDISC,TAU1,TAU2
#endif
C
      SNXT = BIG
      SAFZ1 = PAR(3)-X(3)
      SAFZ2 = PAR(3)+X(3)
C
      A2 = PAR(1)*PAR(1)
      B2 = PAR(2)*PAR(2)
C
      IF(IACT.LT.3)THEN
C
C        -----------------------------------
C        |  Compute safety-distance 'SAFE' |
C        -----------------------------------
C
         X0 = ABS(X(1))
         Y0 = ABS(X(2))
C
         A=PAR(1)
         B=PAR(2)
         X1=X0
         Y1=SQRT((A-X0)*(A+X0))*B/A
         Y2=Y0
         X2=SQRT((B-Y0)*(B+Y0))*A/B
         D1=(X1-X0)**2+(Y1-Y0)**2
         D2=(X2-X0)**2+(Y2-Y0)**2
         DO 1 I=1,8
            IF (B.LT.A) THEN
               X3=.5*(X1+X2)
               Y3=SQRT((A-X3)*(A+X3))*B/A
            ELSE
               Y3=.5*(Y1+Y2)
               X3=SQRT((B-Y3)*(B+Y3))*A/B
            END IF
            IF (D1.LT.D2) THEN
               X2=X3
               Y2=Y3
               D2=(X2-X0)**2+(Y2-Y0)**2
            ELSE
               X1=X3
               Y1=Y3
               D1=(X1-X0)**2+(Y1-Y0)**2
            END IF
    1    CONTINUE
    2    SAFR=SQRT(D1)-1.E-3
*
         SAFE = MIN(SAFZ1,SAFZ2,SAFR)
         IF(IACT.EQ.0)GOTO 99
         IF(IACT.EQ.1.AND.SNEXT.LT.SAFE)GOTO 99
C
      ENDIF
C
C        -----------------------------------
C        |  Compute vector-distance 'SNXT' |
C        -----------------------------------
C
C ....  First check Z
C
      IF(X(6).GT.0.)THEN
         SNXT = SAFZ1/X(6)
      ELSEIF(X(6).LT.0.)THEN
         SNXT = -SAFZ2/X(6)
      ENDIF
C
C ....  Then,if necessary,find the intersection of
C       the given ray(described by array X) whit
C       the cylinder.
C       Ray equation : X'(1-3) = X(1-3) + TAU*X(4-6)
C       Cylinder equation : x**2/a**2 + y**2/b**2 = 1
C       To obtain TAU,solve the quadratic equation
C       Ut**2 + 2Vt + W = 0
C
      SZ = SNXT
      XZ = X(1)+X(4)*SZ
      YZ = X(2)+X(5)*SZ
      IF((XZ*XZ/A2+YZ*YZ/B2).LE.1)GOTO 99
C
      U = X(4)*X(4)*B2+X(5)*X(5)*A2
      V = X(1)*X(4)*B2+X(2)*X(5)*A2
      W = X(1)*X(1)*B2+X(2)*X(2)*A2-A2*B2
      DISCR = V*V-U*W
      IF(DISCR.LT.0.)GOTO 99
      IF(U.EQ.0.)GOTO 99
      SQDISC = SQRT(DISCR)
      TAU1 = (-V+SQDISC)/U
      TAU2 = (-V-SQDISC)/U
C
C ....  Set SNXT to the smallest positive TAU
C
      IF(TAU1.LT.0.)THEN
         IF(TAU2.LT.0.)GOTO 99
         SNXT = TAU2
      ELSE
         SNXT = TAU1
         IF(TAU2.GT.0.0.AND.TAU2.LT.TAU1)THEN
            SNXT = TAU2
         ENDIF
      ENDIF
C
   99 END
