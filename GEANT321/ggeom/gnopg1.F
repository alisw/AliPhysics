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
      SUBROUTINE GNOPG1(X,P,SNXT)
*********************************************************************
***** GNOPG1 ********************************************************
*
*     GNOPG1 ... 20-DEC-1988
*     Version 1.1
*     Rolf Nierhaus
*
*********************************************************************
*
*     Copyright   CERN,   Geneva  1988  -  Copyright  and  any  other
*     appropriate legal protection of  these  computer  programs  and
*     associated  documentation  reserved  in  all  countries  of the
*     world.
*
*********************************************************************
*
*          Subroutine  GNOPG1 is called by GNOPGO for the computation
*     of SNXT, the distance from a point P  along  a  track  T  to  a
*     boundary  surface  of a Geant volume V of shape PGON. The point
*     P is outside the volume V. If the track T has  no  intersection
*     with the volume V, the vector distance SNXT is set to 1.E10.
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
*          GNOPG1 calls GNOPG2 for each section, GNOPG2  call  GNOPG3
*     for  each  sector.  GNOPG4  is  called  to  store  the  surface
*     parameters of  a  sector  in  the  common  block  GCQ1.  GNOPG6
*     computes  the  vector  distance  to  a  convex  volume.  GNOPG7
*     computes the  vector  distance  to  a  plane  surface.  Logical
*     function  GNOPG8  determines  if  a  point  is  inside a convex
*     volume, and logical function GNOPG9 determines if  a  point  is
*     inside a region delimited by a plane surface.
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
*          Logical   function   GNOPG9  returns  TRUE  if  the  point
*     (XP,YP,ZP) is within the positive region of  the  surface  with
*     index   I.   This   is  the  case  if  the  scalar  product  of
*     (XP-XS,YP-YS,ZP-ZS) and (XN,YN,ZN) is positive (or zero).
*
*          GNOPG8 is true if  the  point  (XP,YP,ZP)  is  within  the
*     positive region of all bounding surfaces.
*
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
*     variable  FLAG  TRUE.  GNOPG7  is  called with three parameters
*     I,FLAG and DIST. I is the index of the  surface,  and  DIST  is
*     the vector distance if FLAG is TRUE.
*
*          To  find the vector distance from a point to an elementary
*     volume, all bounding surfaces of the volume are considered.  If
*     the  point  is  in  the positive region of a surface, the track
*     could possibly exit through the surface, but  it  cannot  enter
*     through  it.  For  those surfaces which have the point in their
*     negative region, we  determine  if  the  track  intersects  the
*     surface,  and  what  is the distance to the intersection point.
*     Only the largest of these distances  is  a  candidate  for  the
*     distance  from  the point to the volume. It is however possible
*     that the track misses the elementary volume entirely.  This  we
*     find  out  by applying the function GNOPG8 (Inside volume test)
*     to the coordinates of the intersection point. GNOPG6 is  called
*     with  two parameters: a logical variable FLAG, which is TRUE if
*     the track  intersects  the  volume,  and  DIST,  which  is  the
*     distance if FLAG is TRUE.
*
*          The  coordinates of  the point P and the direction cosines
*     of the  track  T  are  stored  in  the  common  block  GCQ2  by
*     subroutine  GNOPG1. The parameter array P in the call to GNOPG1
*     contains in its first two members the lower  phi-limit  of  the
*     PGON  and  the range in phi. Both angles are in degrees. GNOPG1
*     convertes from degrees to radians and stores the phi-limits  of
*     the  first  sector  of  each  section  in the common block GCQ3
*     together with the number of sectors per section. The number  of
*     sectors  per  section  is  contained in the third member of the
*     parameter array P. The other members of P  have  the  following
*     meaning:  From  P(5)  onwards,  there are groups of three. Each
*     group describes a section boundary. The  first  member  is  the
*     Z-coordinate  of the boundary. The second and the third are the
*     distances from the Z-axis to the inner and outer edges at  that
*     boundary.  If  there is no inner edge, the sector is limited by
*     the Z-axis. In this case the second members of  the  group  are
*     zero.  Groups may be shared by adjacent sections. Otherwise the
*     Z-coordinates of two consecutive  groups  are  the  same.  P(4)
*     contains the number of groups.
*
*          GNOPG1  calls  GNOPG2  with  8 parameters. The first 6 are
*     input parameters, the first  3  refer  to  the  "left"  section
*     boundary,  the next 3 to the "right" section boundary. The last
*     two parameters are output. Logical variable  FLAG  is  TRUE  if
*     the  track  intersects  the  section.  In this case DIST is the
*     distance.
*
*          GNOPG2 calls GNOPG3 with four parameters. The first 2  are
*     input  parameters,  namely  the  phi-limits  of the sector. The
*     last two parameters again are  output:  FLAG  is  TRUE  if  the
*     track   intersects  the  sector.  In  this  case  DIST  is  the
*     distance. The other variables required by GNOPG3 are  the  same
*     for  all  sectors  of the same section and are stored by GNOPG2
*     in the common block GCQ4. GNOPG2  also  stores  in  the  common
*     block GCQ5 the number of boundary surfaces of each sector.
*
*          The  surfaces  orthogonal  to  the Z-axis are the same for
*     all sectors of a section. The surface parameters of  these  two
*     sections  are  stored  by  calls  to  GNOPG4  from  GNOPG2. The
*     surface parameters of the other  three  or  four  surfaces  are
*     stored from GNOPG3.
*
*          GNOPG3  sets  FLAG  TRUE,  if  the  track T intersects the
*     corresponding sector. GNOPG2 finds  the  shortest  distance  to
*     all  intersected  sectors,  and  set  FLAG TRUE, if the track T
*     intersects  the  corresponding  section.   GNOPG1   finds   the
*     shortest  distance  to  all intersected sections. If no section
*     intersects, the track FLAG is set FALSE, and 1.E10 is  returned
*     as SNXT.
*
***** Subroutine GNOPG1 *************************** 04-DEC-1988 *****
#if !defined(CERNLIB_SINGLE)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (F=0.01745329251994330D0)
#endif
#if defined(CERNLIB_SINGLE)
      PARAMETER (F=0.01745329251994330)
#endif
      REAL X(6),P(49),SNXT
      LOGICAL FLAG,FLAG1
      DIMENSION XS(6), YS(6), ZS(6), XN(6), YN(6), ZN(6)
      LOGICAL FLAG1X, GNOP1X, GNOP2X, GNOP4X
      PARAMETER (ONE=1,HALF=ONE/2)
      XP=X(1)
      YP=X(2)
      ZP=X(3)
      XD=X(4)
      YD=X(5)
      ZD=X(6)
      NT=P(3)+.5
      P1=F*P(1)
      P2=F*P(2)/NT
      INDEX=5
      MINDEX=3.*P(4)+1.5
      FLAG=.FALSE.
      DIST=1.E10
   10 Z1=P(INDEX)
      D1N=P(INDEX+1)
      D1X=P(INDEX+2)
      Z2=P(INDEX+3)
      D2N=P(INDEX+4)
      D2X=P(INDEX+5)
C*****  Code Expanded From Routine:  GNOPG2
C*****  Code Expanded From Routine:  GNOPG4
      XS(1) = 0.
      YS(1) = 0.
      ZS(1) = Z1
      XN(1) = 0.
      YN(1) = 0.
      ZN(1) = 1.
C*****  End of Code Expanded From Routine:  GNOPG4
C*****  Code Expanded From Routine:  GNOPG4
      XS(2) = 0.
      YS(2) = 0.
      ZS(2) = Z2
      XN(2) = 0.
      YN(2) = 0.
      ZN(2) = -1.
C*****  End of Code Expanded From Routine:  GNOPG4
      P3 = P1
      P4 = P1 + P2
      Z3 = HALF*(Z1 + Z2)
      D3X = HALF*(D1X + D2X)
      TH1 = ATAN((D2X - D1X)/(Z2 - Z1))
      COSTH1 = COS(TH1)
      SINTH1 = SIN(TH1)
      D3N = HALF*(D1N + D2N)
      ISMAX = 5
      IF (D3N .NE. 0.) THEN
         ISMAX = 6
         TH2 = ATAN((D2N - D1N)/(Z2 - Z1))
         COSTH2 = COS(TH2)
         SINTH2 = SIN(TH2)
      ENDIF
      FLAG1 = .FALSE.
      DIST1 = 1.E10
      DO 60  I = 1, NT
C*****  Code Expanded From Routine:  GNOPG3
C*****  Code Expanded From Routine:  GNOPG4
         XS(3) = 0.
         YS(3) = 0.
         ZS(3) = Z3
         XN(3) = -SIN(P3)
         YN(3) = COS(P3)
         ZN(3) = 0.
C*****  End of Code Expanded From Routine:  GNOPG4
C*****  Code Expanded From Routine:  GNOPG4
         XS(4) = 0.
         YS(4) = 0.
         ZS(4) = Z3
         XN(4) = SIN(P4)
         YN(4) = -COS(P4)
         ZN(4) = 0.
C*****  End of Code Expanded From Routine:  GNOPG4
         P1X = HALF*(P3 + P4)
         COSP = COS(P1X)
         SINP = SIN(P1X)
C*****  Code Expanded From Routine:  GNOPG4
         XS(5) = D3X*COSP
         YS(5) = D3X*SINP
         ZS(5) = Z3
         XN(5) = -COSP*COSTH1
         YN(5) = -SINP*COSTH1
         ZN(5) = SINTH1
C*****  End of Code Expanded From Routine:  GNOPG4
         IF (D3N .NE. 0.) THEN
C*****  Code Expanded From Routine:  GNOPG4
            XS(6) = D3N*COSP
            YS(6) = D3N*SINP
            ZS(6) = Z3
            XN(6) = COSP*COSTH2
            YN(6) = SINP*COSTH2
            ZN(6) = -SINTH2
C*****  End of Code Expanded From Routine:  GNOPG4
         ENDIF
C*****  Code Expanded From Routine:  GNOPG6
         FLAG1X = .FALSE.
         DIST2X = 0.
         DO 20  IS = 1, ISMAX
C*****  Code Expanded From Routine:  GNOPG9
*     TRUE if (XP,YP,ZP) in positive region of surface I
            GNOP1X = 0. .LE. (XP-XS(IS))*XN(IS)+(YP-YS(IS))*YN(IS)+(ZP-
     +         ZS(IS))*ZN(IS)
C*****  End of Code Expanded From Routine:  GNOPG9
C*****  Code Expanded From Routine:  GNOPG7
            IF (.NOT.GNOP1X) THEN
               SPPMSN = (XP - XS(IS))*XN(IS) + (YP - YS(IS))*YN(IS) + (
     +            ZP - ZS(IS))*ZN(IS)
               SPDN = XD*XN(IS) + YD*YN(IS) + ZD*ZN(IS)
               IF (SPDN .EQ. 0.) THEN
                  DIST1X = -.000001
               ELSE
                  DIST1X = -SPPMSN/SPDN
               ENDIF
C*****  End of Code Expanded From Routine:  GNOPG7
               IF ((-1.E-5) .LE. DIST1X) THEN
                  FLAG1X = .TRUE.
                  DIST2X = MAX(DIST1X,DIST2X)
               ENDIF
            ENDIF
   20    CONTINUE
         IF (.NOT.FLAG1X) GO TO 50
         T = DIST2X + .001
         XQ = XP + T*XD
         YQ = YP + T*YD
         ZQ = ZP + T*ZD
C*****  Code Expanded From Routine:  GNOPG8
*     TRUE if (XP,YP,ZP) in volume
         GNOP2X = .FALSE.
         DO 30  IS1X = 1, ISMAX
C*****  Code Expanded From Routine:  GNOPG9
*     TRUE if (XP,YP,ZP) in positive region of surface I
            GNOP4X = 0. .LE. (XQ-XS(IS1X))*XN(IS1X)+(YQ-YS(IS1X))*YN(
     +         IS1X)+(ZQ-ZS(IS1X))*ZN(IS1X)
            IF (.NOT.GNOP4X) GO TO 40
C*****  End of Code Expanded From Routine:  GNOPG9
   30    CONTINUE
         GNOP2X = .TRUE.
   40    CONTINUE
         FLAG1X = GNOP2X
C*****  End of Code Expanded From Routine:  GNOPG8
   50    CONTINUE
C*****  End of Code Expanded From Routine:  GNOPG3
         IF (FLAG1X) THEN
            FLAG1 = .TRUE.
            DIST1 = MIN(DIST2X,DIST1)
         ENDIF
         P3 = P4
         P4 = P4 + P2
   60 CONTINUE
C*****  End of Code Expanded From Routine:  GNOPG2
      IF (FLAG1) THEN
         FLAG=.TRUE.
         IF (DIST1.LT.DIST) DIST=DIST1
      END IF
      INDEX=INDEX+3
      IF (MINDEX.LT.INDEX) THEN
         IF(FLAG) THEN
            SNXT=DIST
         ELSE
            SNXT=1.E10
         ENDIF
      ELSE
         IF (P(INDEX+3).EQ.Z2) INDEX=INDEX+3
         GO TO 10
      ENDIF
      END
