*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:20:21  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/03 29/08/94  11.37.51  by  S.Giani
*-- Author :
*
      SUBROUTINE GDCGVW (VPAR,VMA)
C.
C.    ******************************************************************
C.    *                                                                *
C.    *  Subroutine evaluates the Viewing Matrix from Viewing Angles   *
C.    *                                                                *
C.    *  Input Parameters:                                             *
C.    *                                                                *
C.    *     VPAR(3) : Viewing Angles (VPAR(1)=THETA, VPAR(2)=PHI,      *
C.    *               VPAR(3)=PSI)                                     *
C.    *                                                                *
C.    *  Output Parameters:                                            *
C.    *                                                                *
C.    *     VMA(4,3): Viewing Matrix                                   *
C.    *                                                                *
C.    *    ==>Called by :   GDINIT,GDRAW                               *
C.    *                                                                *
C.    *       Author : J.Salt     *********                            *
C.    *                                                                *
C.    ******************************************************************
C.
#include "geant321/gcdraw.inc"
#include "geant321/gconsp.inc"
*
      DIMENSION  VPAR(3),AMTR(4,4),VMA(4,3)
      DIMENSION  VL(3),VM(3),VN(3)
      COMMON/PROSP/SVN(3)
*.______________________________________
*
*      CALL ISCLIP(1)
      IERR=0
      DO 20 KK=1,4
         DO 10 JJ=1,4
            AMTR(KK,JJ)=0.
   10    CONTINUE
   20 CONTINUE
*
*         Normalise theta, phi and psi angles to [0-360] range
*
      THETA = MOD(ABS(VPAR(1)),180.)
      PHI = MOD(ABS(VPAR(2)),360.)
      PSI = MOD(ABS(VPAR(3)),360.)
*
*         Ensure theta is in the range [0-180]
*
      IF (THETA.GT.180.) THEN
         PHI = PHI + 180.
         THETA = 360. - THETA
      ENDIF
      SINPSI = SIN(PSI * DEGRAD)
      COSPSI = COS(PSI * DEGRAD)
      ST = SIN(THETA * DEGRAD)
      CT = COS(THETA * DEGRAD)
      SP = SIN(PHI * DEGRAD)
      CP = COS(PHI * DEGRAD)
*
*         VN is new nu axis
*
      VN(1) = ST * CP
      SVN(1)=VN(1)
      VN(2) = ST * SP
      SVN(2)=VN(2)
      VN(3) = CT
      SVN(3)=VN(3)
*
      IF (ABS(VN(2)).LT.0.99999) THEN
*
*         Y-axis is default mu axis (view up vector)
*
         VM(1) = 0.
         VM(2) = 1.
         VM(3) = 0.
*
*           Define new lambda axis
*
         VL(1)=VM(2)*VN(3)-VM(3)*VN(2)
         VL(2)=VM(3)*VN(1)-VM(1)*VN(3)
         VL(3)=VM(1)*VN(2)-VM(2)*VN(1)
         VAVL=SQRT(VL(1)**2+VL(2)**2+VL(3)**2)
         VL(1)=VL(1)/VAVL
         VL(2)=VL(2)/VAVL
         VL(3)=VL(3)/VAVL
*
*           Define new mu axis
*
         VM(1)=VN(2)*VL(3)-VN(3)*VL(2)
         VM(2)=VN(3)*VL(1)-VN(1)*VL(3)
         VM(3)=VN(1)*VL(2)-VN(2)*VL(1)
*
      ELSE
*
*           Special case when observer line of sight is along mu:
*           in this case one chooses arbitrarily the vertical axis of the
*           plane of projection as the lambda axis and the horizontal axis
*           as the nu axis
*
         VL(1) = 0.
         VL(2) = 0.
         VL(3) = 1.
         VM(1) = 1.
         VM(2) = 0.
         VM(3) = 0.
      ENDIF
*
*         Get the view up vector by rotating the mu axis
*         PSI degrees in the view plane (= mu-nu plane)
*
      VUPX = COSPSI * VM(1) + SINPSI * VL(1)
      VUPY = COSPSI * VM(2) + SINPSI * VL(2)
      VUPZ = COSPSI * VM(3) + SINPSI * VL(3)
*
*      Values for View Reference Point in NDC (CSw =1)
*
      VRPX = 0.5
      VRPY = 0.5
      VRPZ = 0.5
      VPNX = VN(1)
      VPNY = VN(2)
      VPNZ = VN(3)
*
*           Evaluate view matrix
*
      CALL GDCGEM(VRPX,VRPY,VRPZ,VUPX,VUPY, VUPZ,VPNX, VPNY,VPNZ,IERR,
     +AMTR)
*
      IF (IERR.NE.0)RETURN
*
      VMA(1,1)=AMTR(1,1)*GSCU
      VMA(2,1)=AMTR(1,2)*GSCU
      VMA(3,1)=AMTR(1,3)*GSCU
      VMA(4,1)=GU0
*
      VMA(1,2)=AMTR(2,1)*GSCV
      VMA(2,2)=AMTR(2,2)*GSCV
      VMA(3,2)=AMTR(2,3)*GSCV
      VMA(4,2)=GV0
*
      DO 30  K=1,4
         VMA(K,3)=AMTR(3,K)
   30 CONTINUE
*
      END
