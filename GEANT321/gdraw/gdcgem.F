*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:20:20  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.25  by  S.Giani
*-- Author :
*
      SUBROUTINE GDCGEM(VRPX,VRPY,VRPZ,VUPX,VUPY,VUPZ,VPNX, VPNY,VPNZ,
     +IERR,AMTR)
C.
C.    ******************************************************************
C.    *                                                                *
C.    *                                                                *
C.    *  This subroutine evaluates the Viewing Matrix with the same    *
C.    *  algorithm and with the same input parameters as subroutine    *
C.    *  GEVVM of GKS3D.                                               *
C.    *                                                                *
C.    *  Input Parameters:                                             *
C.    *                                                                *
c.    *       VRPX,VRPY,VRPZ,VUPX,VUPY,VUPZ,VPNX, VPNY,VPNZ            *
C.    *                                                                *
C.    *  Output Parameters:                                            *
C.    *                                                                *
C.    *       IERR,AMTR(4,4)                                           *
C.    *                                                                *
C.    *                                                                *
C.    *       Called by :  GDCGVW                                      *
C.    *                                                                *
C.    *       Author : J.Salt     *********                            *
C.    *                                                                *
C.    ******************************************************************
C.
#include "geant321/gcunit.inc"
*
      DIMENSION VRPT(3),VUPT(3),VPNT(3)
      DIMENSION AMTR(4,4)
*
      IERR=0
      VRPT(1)=VRPX
      VRPT(2)=VRPY
      VRPT(3)=VRPZ
      VUPT(1)=VUPX
      VUPT(2)=VUPY
      VUPT(3)=VUPZ
      VPNT(1)=VPNX
      VPNT(2)=VPNY
      VPNT(3)=VPNZ
*
      VLEN=SQRT(VPNT(1)**2+VPNT(2)**2+VPNT(3)**2)
      AMTR(3,1)=VPNT(1)/VLEN
      AMTR(3,2)=VPNT(2)/VLEN
      AMTR(3,3)=VPNT(3)/VLEN
*
*    First Row (U-axis)  := Unit Vector Perpendicular to VUP and VPN
*
      UX=VUPT(2)*VPNT(3)-VUPT(3)*VPNT(2)
      UY=VUPT(3)*VPNT(1)-VUPT(1)*VPNT(3)
      UZ=VUPT(1)*VPNT(2)-VUPT(2)*VPNT(1)
      ULEN=SQRT(UX**2+UY**2+UZ**2)
*
*     LENGTH=0. <=> VUP and VPN are collinear
*
      IF(ULEN.EQ.0)THEN
         WRITE(CHMAIL,10000)
         CALL GMAIL(0,0)
      ELSE
         AMTR(1,1)=UX/ULEN
         AMTR(1,2)=UY/ULEN
         AMTR(1,3)=UZ/ULEN
      ENDIF
*
*       2ND ROW(V-AXIS)  := UNIT VECTOR PERPENDICULAR TO U- AND N-AXIS
*
      AMTR(2,1)=AMTR(3,2)*AMTR(1,3)-AMTR(1,2)*AMTR(3,3)
      AMTR(2,2)=AMTR(3,3)*AMTR(1,1)-AMTR(1,3)*AMTR(3,1)
      AMTR(2,3)=AMTR(3,1)*AMTR(1,2)-AMTR(1,1)*AMTR(3,2)
*
*      4TH COLUMN : = TRANSLATION FACTORS
*
      AMTR(1,4)=-(VRPT(1)*AMTR(1,1)+VRPT(2)*AMTR(1,2)+ VRPT(3)*
     +AMTR(1,3))
      AMTR(2,4)=-(VRPT(1)*AMTR(2,1)+VRPT(2)*AMTR(2,2)+ VRPT(3)*
     +AMTR(2,3))
      AMTR(3,4)=-(VRPT(1)*AMTR(3,1)+VRPT(2)*AMTR(3,2)+ VRPT(3)*
     +AMTR(3,3))
*
*      4TH ROW
*
      AMTR(4,1)=0.
      AMTR(4,2)=0.
      AMTR(4,3)=0.
      AMTR(4,4)=0.
*
10000 FORMAT(' **** VUP and VPN are collinear  ****')
*
      END
