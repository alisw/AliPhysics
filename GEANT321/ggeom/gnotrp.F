*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:20:53  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.30  by  S.Giani
*-- Author :
      SUBROUTINE GNOTRP(X,P,IACT,SNEXT,SNXT,SAFE)
C.
C.    ******************************************************************
C.    *                                                                *
C.    *     Compute distance up to intersection with 'TRAP' volume.    *
C.    *     from outside point X(1-3) along direction X(4-6)           *
C.    *                                                                *
C.    *     P    (input) : volume parameters                           *
C.    *     IACT (input) : action flag                                 *
C.    *         = 0  Compute SAFE only                                 *
C.    *         = 1  Compute SAFE, and SNXT only if SNEXT.GT.SAFE      *
C.    *         = 2  Compute both SAFE and SNXT                        *
C.    *         = 3  Compute SNXT only                                 *
C.    *     SNEXT (input)    : see IACT = 1                            *
C.    *     SNXT  (output)   : distance to volume boundary             *
C.    *     SAFE  (output)   : shortest distance to any boundary       *
C.    *                                                                *
C.    *     Scalar Distance to Volume.                                 *
C.    *     The scalar distance from a point to a hexahedron can be    *
C.    *     the distance from the point to a surface,                  *
C.    *     the distance from the point to an edge or                  *
C.    *     the distance from the point to a vertex.                   *
C.    *     Here we compute only the distances to the planes of the    *
C.    *     six surfaces and take the maximum.                         *
C.    *     Since the distances to edges or vertices can only be       *
C.    *     larger, this is a first approximation to                   *
C.    *     SAFETY, that is a value which is not larger than the       *
C.    *     distance from the point to the volume.                     *
C.    *                                                                *
C.    *    ==>Called by : GNEXT, GTNEXT                                *
C.    *         Author  R.Nierhaus  *********                          *
C.    *                                                                *
C.    ******************************************************************
*
#include "geant321/gconsp.inc"
*
      DIMENSION X(6),P(35)
*
      SAFE=-BIG
      SNXT=BIG
      DO 1 I=12,32,4
         T1=P(I)*X(1)+P(I+1)*X(2)+P(I+2)*X(3)+P(I+3)
         IF (SAFE.LT.T1) SAFE=T1
    1 CONTINUE
      IF (IACT.EQ.0.OR.(IACT.EQ.1.AND.SNEXT.LE.SAFE)) RETURN
*
*     Vector Distance to Volume.
*     Volume is hexahedron.
*     (X(1),X(2),X(3)) is outside.
*     P(1),P(2) ... P(11) are the standard Geant execution time
*     parameters for a shape 4 (TRAP) volume.
*     P(12) ... P(35) contain boundary surface data.
*     P(12),P(13),P(14),P(15) are the coefficients of the normalized
*     implicit plane equation for the first boundary surface.
*     P(16),P(17),P(18),P(19) for the second boundary surface.
*     P(32),P(33),P(34),P(35) for the last (sixth) boundary surface.
*     The first three coefficients are the components of a unit
*     vector pointing away from the volume.
*     T1 is the signed orthogonal distance of the point
*     (X(1),X(2),X(3)) and a boundary plane.
*     If this distance is negative, the track cannot enter the
*     volume through the corresponding face.
*     T2 is the cosine of the angle between the plane normal and the
*     track directions (X(4),X(5),X(6)).
*     T3 is the vector distance to the corresponding surface.
*
      SNXT1=-BIG
      DO 2 I=12,32,4
         T1=P(I)*X(1)+P(I+1)*X(2)+P(I+2)*X(3)+P(I+3)
         IF (0.LT.T1) THEN
            T2=P(I)*X(4)+P(I+1)*X(5)+P(I+2)*X(6)
            T3=-.0000001
            IF (T2.NE.0.) T3=-T1/T2
            IF (0..LE.T3) THEN
               IF (SNXT1.LT.T3) SNXT1=T3
            END IF
         END IF
    2 CONTINUE
      IF (SNXT1.EQ.-BIG) RETURN
*
      T=1.00001*SNXT1+.00001
*
*     (XQ,YQ,ZQ) is the intersection point of the track with
*     a boundary surface plane.
*     P(I)*XQ+P(I+1)*YQ+P(I+2)*ZQ+P(I+3) is the signed distance of
*     the intersection point with a boundary surface plane.
*     If this distance is positive for any of the six surfaces,
*     the intersection is not with the volume.
*     A small correction is applied which moves the point slightly
*     to the interior of the volume to protect against rounding
*     errors.
*
      XQ=X(1)+T*X(4)
      YQ=X(2)+T*X(5)
      ZQ=X(3)+T*X(6)
      DO 3 I=12,32,4
         IF (0.LT.P(I)*XQ+P(I+1)*YQ+P(I+2)*ZQ+P(I+3)) RETURN
    3 CONTINUE
      SNXT=SNXT1
*
      END
