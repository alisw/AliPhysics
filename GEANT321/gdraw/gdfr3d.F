*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:20:23  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.26  by  S.Giani
*-- Author :
      SUBROUTINE GDFR3D(X,NPOINT,U,V)
C.
C.    ******************************************************************
C.    *                                                                *
C.    *       Routine to convert the N 3D points given by array X      *
C.    *       to the 2D points in U and V                              *
C.    *                                                                *
C.    *       If IGMR=1 then 3D points are not projected onto          *
C.    *       U-V plane but used to build the APOLLO-GMR 3D structure  *
C.    *                                                                *
C.    *       ITR3D  = 0 for standard projection;                      *
C.    *             <> 0 for rotation + projection i.e. R-Z projection *
C.    *                 with sign of R located in ITRSGN [ITR3D]       *
C.    *                                                                *
C.    *                                                                *
C.    *       Steps from 3D to 2D when NPOINT>0 :                      *
C.    *                                                                *
C.    *       a) for volumes (IOBJ=1) :                                *
C.    *          transformation of the 3D point (x',y',z') in DRS      *
C.    *          (Daughter Reference System) into the 3D point (x,y,z) *
C.    *          in MARS (MAster Reference System);                    *
C.    *                                                                *
C.    *          for tracks/hits (IOBJ<>1) after GDRVOL with NRS<>0 :  *
C.    *          transformation of the 3D point (x',y',z') in MARS     *
C.    *          into the 3D point (x,y,z) in DRS;                     *
C.    *                                                                *
C.    *          for tracks/hits (IOBJ<>1) after GDRAW/GDRAWC :   *
C.    *          nothing                                               *
C.    *                                                                *
C.    *       b) transformation of the 3D point (x,y,z) in MARS/DRS    *
C.    *          into the 3D point (xx,yy,zz) in PRS (Projection       *
C.    *          Reference System); PRS has its z axis along           *
C.    *          the observer's line of sight (going to into the eye)  *
C.    *                                                                *
C.    *       c) transformation of the 3D point (xx,yy,zz) in PRS      *
C.    *          into the 2D point (u,v) in the plane normal to the    *
C.    *          observer's line of sight; at that stage if ITR3D=0    *
C.    *          (standard projection) then u=xx and v=yy,             *
C.    *          but if ITR3D <> 0 (rotation + projection i.e. R-Z)    *
C.    *          then the 3D point (xx,yy,zz) in PRS is rotated along  *
C.    *          xx axis until it is onto the xx-yy plane and finally  *
C.    *          u=xx (unchanged by the rotation) and v=yy             *
C.    *                                                                *
C.    *    ==>Called by : <USER>, GDAHIT, GDARC, GDAXIS, GDCIRC,       *
C.    *                   GDCIRR, GDCIRY, GDCUT, GDCXYZ, GDLINE,       *
C.    *                   GDPART, GDRECT, GDSURF, GDXYZ                *
C.    *       Authors : R.Brun, P.Zanarini   *********                 *
C.    *                                                                *
C.    ******************************************************************
C.
#include "geant321/gcvolu.inc"
#include "geant321/gcdraw.inc"
      DIMENSION XC(3),XYZ(3),VL(3),VM(3),VN(3),X(3,1),U(1),V(1)
      DIMENSION XP(3,50)
      SAVE VL,VM,VN
C
C.
C.    ------------------------------------------------------------------
C.
      CALL UCTOH('PERS',IPERS,4,4)
C
      IF(NPOINT.EQ.0)GO TO 100
      N=NPOINT
      IF(NPOINT.LT.0)N=-NPOINT
      IF(NGVIEW.NE.0) GO TO 40
C
C             First call compute the rotation matrix
C
      PH     = ABS(MOD(GPHI,360.))
      THET   = ABS(MOD(GTHETA,360.))
      IF(THET.LE.180.)GO TO 10
      PH     = PH + 180.
      THET   = 360. - THET
C
   10 ST     = SIN(THET * 0.017453)
      CT     = COS(THET * 0.017453)
      SP     = SIN(PH * 0.017453)
      CP     = COS(PH * 0.017453)
C
C             VN is new nu axis
C
      VN(1)  = ST * CP
      VN(2)  = ST * SP
      VN(3)  = CT
C
      IF(ABS(VN(2)).GT.0.99999)GO TO 20
C
      VM(1)  = 0.
      VM(2)  = 1.
      VM(3)  = 0.
C
C             Define new lambda axis
C
      CALL CROSS(VM,VN,VL)
      CALL VUNIT(VL,VL,3)
C
C             Define new mu axis
C
      CALL CROSS(VN,VL,VM)
      GO TO 30
C
C             Special case when observer line of sight is along mu:
C             in this case one chooses arbitrarily the vertical axis of
C             plane of projection as the lambda axis and the horizontal
C             as the nu axis
C
   20 VL(1)  = 0.
      VL(2)  = 0.
      VL(3)  = 1.
      VM(1)  = 1.
      VM(2)  = 0.
      VM(3)  = 0.
   30 CONTINUE
C
      NGVIEW=1
C
   40 CONTINUE
C
C             Begin of a call with viewing tranformation unchanged (NGVI
C
      DO 90  I=1,N
C
         IF (NPOINT.LT.0) THEN
C
C             NPOINT < 0 : X is in MARS
C
            IF (IGMR.EQ.1) THEN
               DO 50  J=1,3
   50          XP(J,I)=X(J,I)
               GO TO 90
            ENDIF
C
C             Transform X (MARS) in XYZ (PRS)
C             i.e. project onto U,V,W
C
            XYZ(1)=X(1,I)*VL(1)+X(2,I)*VL(2)+X(3,I)*VL(3)
            XYZ(2)=X(1,I)*VM(1)+X(2,I)*VM(2)+X(3,I)*VM(3)
            XYZ(3)=X(1,I)*VN(1)+X(2,I)*VN(2)+X(3,I)*VN(3)
C
         ELSE
C
            IF (IOBJ.NE.1) THEN
C
C             NPOINT > 0 and IOBJ <> 1 : X belongs to a track or hit or
C             so it must be transformed from MARS to last DRS used by GD
C             (otherwise GTRAN0 and GRMAT0 are the unitary transformatio
C             and X is still expressed in MARS)
C
               CALL GTRNSF(X(1,I),GTRAN0,GRMAT0,XC)
C
            ELSE
C
C             NPOINT > 0 and IOBJ = 1 : X belongs to a volume,
C             so it must be transformed from DRS to MARS
C
               CALL GINROT(X(1,I),GRMAT(1,NLEVEL),XC)
               DO 60  J=1,3
   60          XC(J)=XC(J)+GTRAN(J,NLEVEL)
C
            ENDIF
C
            IF (IGMR.EQ.1) THEN
               DO 70  J=1,3
   70          XP(J,I)=XC(J)
               GO TO 90
            ENDIF
C
C             Transform XC (MARS or DRS) in XYZ (PRS)
C             i.e. project onto U,V,W
C
            XYZ(1)=XC(1)*VL(1)+XC(2)*VL(2)+XC(3)*VL(3)
            XYZ(2)=XC(1)*VM(1)+XC(2)*VM(2)+XC(3)*VM(3)
            XYZ(3)=XC(1)*VN(1)+XC(2)*VN(2)+XC(3)*VN(3)
C
         ENDIF
C
C             Scale, shift and rotate in in U and V
C
   80    CONTINUE
C
C             R-Z projection ?
C
         IF (ITR3D.NE.0) CALL GDTHRZ(XYZ)
C
C             Perspective projection ?
C
         IF (IPRJ.EQ.IPERS) THEN
            IF (XYZ(3).GE.DPERS) XYZ(3)=DPERS-0.0001
            XYZ(1)=XYZ(1)*DPERS/(DPERS-XYZ(3))
            XYZ(2)=XYZ(2)*DPERS/(DPERS-XYZ(3))
         ENDIF
C
         UU=XYZ(1)*GSCU
         VV=XYZ(2)*GSCV
         U(I)=GU0+UU*COSPSI-VV*SINPSI
         V(I)=GV0+UU*SINPSI+VV*COSPSI
C
C             Zoom processing
C
         U(I)=U(I)*GZUA+GZUB+GZUC
         V(I)=V(I)*GZVA+GZVB+GZVC
C
   90 CONTINUE
C
      IF (IGMR.EQ.1) CALL GM3POL(XP,N)
C
  100 CONTINUE
      END
