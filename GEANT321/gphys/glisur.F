*PATCH,GEXAM1.
*CMZ :  3.21/02 29/03/94  15.41.35  by  S.Giani
*DECK,GLISUR
*CMZ :          12/09/95  11.22.49  by  S.Ravndal
*CMZ :  3.21/02 03/07/94  19.08.44  by  S.Giani
*-- Author :
      SUBROUTINE GLISUR(X0,X1,MEDI0,MEDI1,U,PDOTU,IERR)
C.
C.    ******************************************************************
C.    *                                                                *
C.    *  This routine simulates the surface profile of a boundary      *
C.    *  between two media as seen by an approaching particle with     *
C.    *  coordinates and direction given by X0. The surface is         *
C.    *  identified by the arguments MEDI0 and MEDI1 which are the     *
C.    *  media indices of the region in which the track is presently   *
C.    *  and the one which it approaches, respectively. The input      *
C.    *  vector X1 contains the coordinates of a point on the other    *
C.    *  side of the boundary from X0, and lying within medium MEDI1.  *
C.    *  The result is a unit vector U normal to the surface of        *
C.    *  reflection at X0 and pointing into the medium from which the  *
C.    *  track is approaching. The quality of the surface finish is    *
C.    *  described by the parameter POLISH, which varies from 0 to 1.  *
C.    *  POLISH=0 means maximum roughness, with effective plane of     *
C.    *  reflection distributed as cos(alpha), where alpha is the      *
C.    *  angle between the unit normal to the effective plane of       *
C.    *  reflection and the normal to the nominal medium boundary at   *
C.    *  X0. POLISH=1 means perfect smoothness. In between the surface *
C.    *  is modeled by a bell-shaped distribution in alpha, with       *
C.    *  limits defined by                                             *
C.    *              sin(alpha) = +/- (1-POLISH)                       *
C.    *  At the interface between two media, the surface POLISH        *
C.    *  parameter is taken from a user routine GUPLSH(MEDI0,MEDI1).   *
C.    *  The indices MEDI0 and MEDI1 refer to the media declared by    *
C.    *  the user. If GGPERP returns an error, there was a precision   *
C.    *  problem with the tracking, and point X0 is on the same side   *
C.    *  of the surface as X1. In this case (or any other error        *
C.    *  condition from GGPERP) a return is made with IERR ^=0. If     *
C.    *  IERR=0 on return then U contains a valid unit vector.         *
C.    *                                                                *
C.    *    ==>Called by : <USER>, UGINIT                               *
C.    *       Author    F.Carminati, R.Jones ***********               *
C.    *                                                                *
C.    ******************************************************************
C.
#include "geant321/gcvolu.inc"
#include "geant321/gcvol1.inc"
C
      REAL X0(6),X1(3),U(3)
      REAL D(3),RNDM(3)

*
* *** Decide which volume defines the surface
      IF (NLEVEL.GE.NLEVL1) THEN
         LVOL = 1
      ELSE
         LVOL = -1
         DO 10 I=2,NLEVEL
            IF ((NAMES1(I).NE.NAMES(I)).OR. (NUMBR1(I).NE.NUMBER(I)))
     +      LVOL = 1
   10    CONTINUE
      END IF
*
* *** Find the nominal unit normal to the surface
      IF (LVOL.EQ.1) THEN
C
C Instead of X0, now X1 is used by GGPERP, in order to calculate
C the right perpendicular even near to a corner,
C see the modifactions in GTCKOV concerning VBOU
C
         CALL GGPERP(X1,U,IERR)
      ELSE
         LVOL = NLEVEL
         NLEVEL = NLEVL1
         CALL GGPERP(X1,U,IERR)
         U(1) = -U(1)
         U(2) = -U(2)
         U(3) = -U(3)
         NLEVEL = LVOL
      END IF
      IF (IERR.NE.0) GO TO 999
      PDOTU = X0(4)*U(1) + X0(5)*U(2) + X0(6)*U(3)
      POLISH=GUPLSH(MEDI0,MEDI1)
      IF (POLISH.LT.1.) THEN
*
* *** Generate distortion vector D within a uniform sphere
         DIAM = 2.*(1.-POLISH)
         DIA2 = DIAM**2
         IROUND = 0
   20    IF(IROUND.GT.5) GO TO 999
         IROUND = IROUND+1
         CALL GRNDM(RNDM,3)
         D(1) = RNDM(1)-0.5
         D(2) = RNDM(2)-0.5
         D(3) = RNDM(3)-0.5
         D2 = D(1)**2+D(2)**2+D(3)**2
         IF (D2.GT.0.25) GO TO 20
         D(1) = D(1)*DIAM
         D(2) = D(2)*DIAM
         D(3) = D(3)*DIAM
         D2 = D2*DIA2
*
* *** Insure that V=U+D will cause reflection away from surface
         PDOTD = X0(4)*D(1) + X0(5)*D(2) + X0(6)*D(3)
         PDOTV = PDOTU+PDOTD
         UDOTD = U(1)*D(1) + U(2)*D(2) + U(3)*D(3)
         V2 = 1+D2+2*UDOTD
         IF (PDOTD*V2+PDOTV*(1.-D2).GT.0.) GO TO 20
*
* *** Normalize V and call it U
         VABS = 1./SQRT(V2)
         U(1) = (U(1)+D(1))*VABS
         U(2) = (U(2)+D(2))*VABS
         U(3) = (U(3)+D(3))*VABS
         PDOTU = PDOTV*VABS
      END IF
  999 END
