*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:20:57  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.31  by  S.Giani
*-- Author :
*
      SUBROUTINE GVLONG (DIR, C, NC, TAU, NTAU)
************************************************************************
*                                                                      *
* GVLONG calculates the distance to the surface along a      VP 920125 *
*        given DIRection                                               *
*                                                                      *
*        SHOULD BE called ONLY AFTER call GVSAFE.                      *
*        due to some results of GVSAFE are used.                       *
*        Just by this reason starting XYZ point                        *
*        is omitted                                                    *
*        =======================================                       *
*                                                                      *
* Input  : DIR -   direction                                           *
*          C   -   coeffs of surface (C(1) not used)                   *
*          NC  -   number of coefs                                     *
*                                                                      *
* Output : TAU -   TAU(1:2) distances to surface in increasing order   *
*                  in the units of the length of DIR                   *
*          NTAU -  number of distances (solutions)  (0,1,2)            *
*                                                                      *
************************************************************************
      COMMON /SLATE/ SLATE (40)
      INTEGER        ISLATE(40)
      EQUIVALENCE    (SLATE,ISLATE)
      REAL           C(*), DIR(3), TAU(2), CBA(3)
      LOGICAL        XLINE
*-----------------------------------------------------------------------
*
*     SLATE(10) is the value of S(X) = equation of the surface with
*     the point coordinates
*     SLATE(11),SLATE(12),SLATE(13) is the gradient to the surface
*
      NTAU = 0
      CBA(1) = SLATE(10)
      CBA(2) = SLATE(11)*DIR(1)+SLATE(12)*DIR(2)+SLATE(13)*DIR(3)
*
      CBA(3) = 0.
      XLINE  = .TRUE.
*
**              in case of a plane
      IF (NC.EQ.4)              GO TO 200
*
*     in case it is a simplified surface X=-C0,Y=-C0,Z=-C0, or X*2+Y*2=-C0
*
      IF (NC.EQ.2)              THEN
        IAX  = C(2)
        IF (IAX.NE.4)           GO TO 200
        CBA(3) = DIR(1)**2+DIR(2)**2
      ELSE
*
**              in case it is a surface with 7 or 10 coefficients
        CBA(3) = DIR(1)**2*C(5)+DIR(2)**2*C(6)+DIR(3)**2*C(7)
*
**              in case it is a surface with 10 coefficients
        IF (NC.EQ.10)           THEN
          CBA(3) = CBA(3) + C(08)*DIR(1)*DIR(2)
          CBA(3) = CBA(3) + C(09)*DIR(2)*DIR(3)
          CBA(3) = CBA(3) + C(10)*DIR(3)*DIR(1)
        ENDIF
      ENDIF
      XLINE = .FALSE.
*
  200 IF (XLINE)                THEN
*
*     in case it is a plane (4 coefficients) or a simplified surface
*
        IF (1.E+6*ABS(CBA(2)).LE.ABS(CBA(1)))   GO TO 999
        NTAU = 1
        TAU(1) = - CBA(1)/CBA(2)
      ELSE
*
**              in case it is a surface with 7 or 10 coefficients
        CALL GVPSQR (CBA, TAU, NTAU)
      ENDIF
*
  999 NTAU = MAX(0,NTAU)
      RETURN
      END
