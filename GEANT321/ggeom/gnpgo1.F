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
*********************************************************************
***** GNPGO1 ********************************************************
*
*     GNPGO1 ... 15-AUG-1991
*     Version 1.1
*     Rolf Nierhaus
*
*********************************************************************
*
*     Copyright   CERN,   Geneva  1991  -  Copyright  and  any  other
*     appropriate legal protection of  these  computer  programs  and
*     associated  documentation  reserved  in  all  countries  of the
*     world.
*
*********************************************************************
*
*          Subroutine  GNPGO1 is called by GNPGON for the computation
*     of SNXT, the distance from a point P  along  a  track  T  to  a
*     boundary  surface  of a Geant volume V of shape PGON. The point
*     P is inside the volume V.
*
*          V  is  generally  a composite volume consisting of several
*     sections. The sections have  boundary  surfaces  orthogonal  to
*     the   Z-axis.   Each  section  consists  generally  of  several
*     sectors. Each sector is an  "elementary"  convex  volume.  This
*     package  assumes it is either a hexahedron or a pentahedron. If
*     it is a pentahedron, it has 6 vertices, of  which  two  are  on
*     the  Z-axis.  All  sectors  of  the same section are congruent.
*     Each section has the same number of sectors.
*
*          We  describe each surface by 6 parameters: the first three
*     are   the   coordinates   of   a   point   on    the    surface
*     XS(I),YS(I),ZS(I),  the  other  three are the components of the
*     normal vector of the surface XN(I),YN(I),ZN(I). I is the  index
*     of  the surface. We consider only one sector at a time, and the
*     number of boundary  surfaces  is  never  larger  then  6.  Each
*     surface  divides  the  space  into  two  regions:  the positive
*     region and the negative region. We choose the direction of  the
*     normal  vectors  of the boundary surfaces such that the bounded
*     volume is within the positive region of each surface, that  is,
*     the normal vector is pointing to the inside of the volume.
*
***** Subroutine GNPGO1 *************************** 15-AUG-1991 *****
      SUBROUTINE GNPGO1(X,P,SNXT)
#if !defined(CERNLIB_SINGLE)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
#endif
      REAL X(6),P(49),SNXT
#if defined(CERNLIB_SINGLE)
      PARAMETER (F=0.01745329251994330,TPI=6.283185307179586)
#endif
#if !defined(CERNLIB_SINGLE)
      PARAMETER (F=0.01745329251994330D0,TPI=6.283185307179586D0)
#endif
      PARAMETER (ONE=1,HALF=ONE/2,ZERO=0.)
      DIMENSION XS(6), YS(6), ZS(6), XN(6), YN(6), ZN(6)
      LOGICAL FLAG, FLG
      XP=X(1)
      YP=X(2)
      ZP=X(3)
      XD=X(4)
      YD=X(5)
      ZD=X(6)
      IMAX=P(4)-.5
*     IMAX   ->   number of Z-sections
      JMAX=P(3)+.5
*     JMAX   ->   number of Phi-sectors
      SNXT=0.
   10 CONTINUE
*     Find current elementary volume
      IF (ZP.LE.P(5)) RETURN
*     Current point (XP,YP,ZP) is below first section.
      DO 20 II=1,IMAX
         IF (ZP.LT.P(5+3*II)) THEN
            I=II
            GO TO 30
         ENDIF
   20 CONTINUE
      RETURN
*     Current point (XP,YP,ZP) is above last section.
   30 CONTINUE
      IF (XP.EQ.0..AND.YP.EQ.0.) XP=1.E-20
      PHI=ATAN2(YP,XP)
      IF (PHI.LT.0.) PHI=PHI+TPI
      P1=F*P(1)
      PHI1=PHI-P1
      IF (PHI1.LT.0.) PHI1=PHI1+TPI
      IF (PHI1.GE.TPI) PHI1=PHI1-TPI
      IF (JMAX.EQ.1) THEN
          IF (ABS(PHI1-TPI).LT.1D-7) PHI1=0.
      ENDIF
      J=PHI1*P(3)/(F*P(2))+ONE
      IF (P(2).EQ.360.) THEN
         IF (J.LT.1) THEN
            J=J+JMAX
         ELSEIF (JMAX.LT.J) THEN
            J=J-JMAX
         END IF
      END IF
      IF (JMAX.LT.J.OR.J.LT.1) RETURN
*     Current point is outside Phi-range.
C*****  Code Expanded From Routine:  GNPGO2
*     GNPGO2 finds the vector distance to the boundary surface
*     of the current elementary volume.
*     I is Z-section, J is Phi-sector.
*     GNPGO2 calls GNPGO4 five or six times for the storage of
*     the surface coefficients of its boundary surfaces.
*
      INDEX = 2 + 3*I
      Z1 = P(INDEX)
      D1N = P(INDEX+1)
      D1X = P(INDEX+2)
      Z2 = P(INDEX+3)
      D2N = P(INDEX+4)
      D2X = P(INDEX+5)
      ZM = HALF*(Z1 + Z2)
      P11X = F*(P(1)+(J-1)*P(2)/JMAX)
      P2 = F*(P(1)+J*P(2)/JMAX)
      PP = HALF*(P11X + P2)
      COSP = COS(PP)
      SINP = SIN(PP)
      DMX = HALF*(D1X + D2X)
      DMN = HALF*(D1N + D2N)
      THX = ATAN((D2X - D1X)/(Z2 - Z1))
      COSTHX = COS(THX)
      SINTHX = SIN(THX)
      XNN = -COSP*COSTHX
      YNN = -SINP*COSTHX
C*****  Code Expanded From Routine:  GNPGO4
*     Store surface coefficients
      XS(5) = DMX*COSP
      YS(5) = DMX*SINP
      ZS(5) = ZM
      XN(5) = XNN
      YN(5) = YNN
      ZN(5) = SINTHX
C*****  End of Code Expanded From Routine:  GNPGO4
C*****  Code Expanded From Routine:  GNPGO9
*          Logical   function   GNPGO9  returns  TRUE  if  the  point
*     (XP,YP,ZP) is within the positive region of  the  surface  with
*     index   I.   This   is  the  case  if  the  scalar  product  of
*     (XP-XS,YP-YS,ZP-ZS) and (XN,YN,ZN) is positive (or zero).
      RESULT=(XP-XS(5))*XN(5)+(YP-YS(5))*YN(5)+(ZP-ZS(5))*ZN(5)
      FLG = 0. .LE. RESULT
      IF (.NOT.FLG) GO TO 50
      ISMAX = 5
      IF (DMN .NE. 0.) THEN
         ISMAX = 6
         THN = ATAN((D2N - D1N)/(Z2 - Z1))
         COSTHN = COS(THN)
         SINTHN = SIN(THN)
         XNN = COSP*COSTHN
         YNN = SINP*COSTHN
C*****  Code Expanded From Routine:  GNPGO4
         XS(6) = DMN*COSP
         YS(6) = DMN*SINP
         ZS(6) = ZM
         XN(6) = XNN
         YN(6) = YNN
         ZN(6) = -SINTHN
C*****  End of Code Expanded From Routine:  GNPGO4
C*****  Code Expanded From Routine:  GNPGO9
         RESULT=(XP-XS(6))*XN(6)+(YP-YS(6))*YN(6)+(ZP-ZS(6))*ZN(6)
         FLG = 0. .LE. RESULT
         IF (.NOT.FLG) GO TO 50
      ENDIF
C*****  Code Expanded From Routine:  GNPGO4
      XS(1) = ZERO
      YS(1) = ZERO
      ZS(1) = Z1
      XN(1) = ZERO
      YN(1) = ZERO
      ZN(1) = ONE
C*****  End of Code Expanded From Routine:  GNPGO4
C*****  Code Expanded From Routine:  GNPGO4
      XS(2) = ZERO
      YS(2) = ZERO
      ZS(2) = Z2
      XN(2) = ZERO
      YN(2) = ZERO
      ZN(2) = -ONE
C*****  End of Code Expanded From Routine:  GNPGO4
C*****  Code Expanded From Routine:  GNPGO4
      XS(3) = ZERO
      YS(3) = ZERO
      ZS(3) = ZM
      XN(3) = -SIN(P11X)
      YN(3) = COS(P11X)
      ZN(3) = ZERO
C*****  End of Code Expanded From Routine:  GNPGO4
C*****  Code Expanded From Routine:  GNPGO4
      XS(4) = ZERO
      YS(4) = ZERO
      ZS(4) = ZM
      XN(4) = SIN(P2)
      YN(4) = -COS(P2)
      ZN(4) = ZERO
C*****  End of Code Expanded From Routine:  GNPGO4
C*****  Code Expanded From Routine:  GNPGO5
*     Vector distance to volume boundary
      SNXT1 = 1.E10
      DO 40  IS = 1, ISMAX
C*****  Code Expanded From Routine:  GNPGO7
*          To  find  the  distance  from  a  point (XP,YP,ZP) along a
*     track  with  direction  cosines   (XD,YD,ZD)   to   a   surface
*     (XS,YS,ZS)(XN,YN,ZN),  we  compute  first the scalar product of
*     the  vector  (XS-XP,YS-YP,ZS-ZP)   with   the   normal   vector
*     (XN,YN,ZN),  then  the scalar product of the vectors (XD,YD,ZD)
*     and (XN,YN,ZN).  The  first  scalar  product  is  the  shortest
*     distance  from  the  point  to  the  plane,  the  second scalar
*     product is the cosine of the angle between the  track  and  the
*     plane  normal.  The  quotient  is  the vector distance. If this
*     vector distance is  positive  (or  zero)  we  set  the  logical
*     variable  FLAG  TRUE.  GNPGO7  is  called with three parameters
*     I,FLAG and DIST. I is the index of the  surface,  and  DIST  is
*     the vector distance if FLAG is TRUE.
         SPPMSN = (XP - XS(IS))*XN(IS) + (YP - YS(IS))*YN(IS) + (ZP -
     +   ZS (IS))*ZN(IS)
         SPDN = XD*XN(IS) + YD*YN(IS) + ZD*ZN(IS)
         IF (SPDN .EQ. 0.) THEN
            DIST1 = 0.
         ELSE
            DIST1 = -(SPPMSN + .0001)/SPDN
         ENDIF
         FLAG = 0. .LT. DIST1
C*****  End of Code Expanded From Routine:  GNPGO7
         IF (FLAG) SNXT1 = MIN(DIST1,SNXT1)
   40 CONTINUE
   50 CONTINUE
C*****  End of Code Expanded From Routine:  GNPGO2
      IF (FLG) THEN
         SNXT=SNXT+SNXT1
         XP=XP+SNXT1*XD
         YP=YP+SNXT1*YD
         ZP=ZP+SNXT1*ZD
*     The current point (XP,YP,ZP) is propagated along the track
*     to the boundary of the current elementary volume.
         GO TO 10
      ENDIF
      END
