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
      SUBROUTINE GVGRAD (XYZ, C, NC, GRA)
************************************************************************
*                                                                      *
* GVGRAD calculates the gradient vector of a surface         VP 880314 *
*                                                                      *
* Input  : XYZ     coordinates of the point                            *
*          C(1)       number of non constant coefficients of the       *
*                        surface                                       *
*          C(2),C(3),... non constant coefficients of the surface      *
*          NC      total number of coefficients of the surface         *
*                                                                      *
* Output : GRA     gradient vector of the surface                      *
*                                                                      *
************************************************************************
      REAL           XYZ(3) , C(*), GRA(3)
*-----------------------------------------------------------------------
*
*     case with simplified surface X=C0, Y=C0, Z=C0, X*2+Y*2=C0
*     (happens only when initialisation is done)
*
      IF (NC.EQ.2)              THEN
        IAX = C(2)
        IF (IAX.LE.3)           THEN
          GRA(1) = 0
          GRA(2) = 0
          GRA(3) = 0
          GRA(IAX) = 1
        ELSE
          GRA(1)   = 2.*XYZ(1)
          GRA(2)   = 2.*XYZ(2)
          GRA(3)   = 0.
        ENDIF
        GO TO 999
      ENDIF
*
*     case with surfaces with 4, 7 or 10 coefficients (normal case)
*
      GRA(1) = C(2)
      GRA(2) = C(3)
      GRA(3) = C(4)
      IF (NC.EQ.4)              GO TO 999
*
*     case with surfaces with 7 or 10 coefficients
*
      DO 100 I = 1,3
  100   GRA(I) = GRA(I) + 2.*C(I+4)*XYZ(I)
      IF (NC.NE.10)             GO TO 999
*
*     case with surfaces with 10 coefficients
*
      GRA(1)   = GRA(1) + C(8)*XYZ(2)+C(10)*XYZ(3)
      GRA(2)   = GRA(2) + C(8)*XYZ(1)+C( 9)*XYZ(3)
      GRA(3)   = GRA(3) + C(9)*XYZ(2)+C(10)*XYZ(1)
 
  999 RETURN
      END
