*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:20:19  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.25  by  S.Giani
*-- Author :
      SUBROUTINE GD3D3D(XIN,NPOINT,XOUT,INV)
C.
C.    ******************************************************************
C.    *                                                                *
C.    *       Routine to transform the 3D points given by array XIN    *
C.    *       into the 3D points in XOUT:                              *
C.    *       if INV=0 the Z-axis for XOUT is the line of sight for XIN*
C.    *       if INV=1 the Z-axis for XIN is the line of sight for XOUT*
C.    *                                                                *
C.    *    ==>Called by : <USER>, GDCUT                                *
C.    *       Author : P.Zanarini   *********                          *
C.    *                                                                *
C.    ******************************************************************
C.
#include "geant321/gcvolu.inc"
#include "geant321/gcdraw.inc"
      DIMENSION XC(3),XOUT(3,1),VL(3),VM(3),VN(3),XIN(3,1)
      SAVE VL,VM,VN
C.
C.    ------------------------------------------------------------------
C.
      IF(NPOINT.EQ.0)GO TO 999
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
C             Begin of a normal call (i.e. with NGVIEW=1)
C
      DO 70  I=1,N
C
C             Loop over the points
C
         IF (NPOINT.LT.0) THEN
C
C             NPOINT < 0 : XIN is in general reference system
C
            IF (INV.EQ.1) THEN
               XC(1)=VL(1)*XIN(1,I)+VM(1)*XIN(2,I)+VN(1)*XIN(3,I)
               XC(2)=VL(2)*XIN(1,I)+VM(2)*XIN(2,I)+VN(2)*XIN(3,I)
               XC(3)=VL(3)*XIN(1,I)+VM(3)*XIN(2,I)+VN(3)*XIN(3,I)
            ELSE
               XC(1)=XIN(1,I)*VL(1)+XIN(2,I)*VL(2)+XIN(3,I)*VL(3)
               XC(2)=XIN(1,I)*VM(1)+XIN(2,I)*VM(2)+XIN(3,I)*VM(3)
               XC(3)=XIN(1,I)*VN(1)+XIN(2,I)*VN(2)+XIN(3,I)*VN(3)
            ENDIF
C
            DO 50  J=1,3
               XOUT(J,I)=XC(J)
   50       CONTINUE
C
         ELSE
C
C             NPOINT > 0 : XIN is in volume reference system
C
            CALL GINROT(XIN(1,I),GRMAT(1,NLEVEL),XC)
            DO 60  J=1,3
               XOUT(J,I)=XC(J)+GTRAN(J,NLEVEL)
   60       CONTINUE
C
         ENDIF
C
   70 CONTINUE
C
  999 RETURN
      END
