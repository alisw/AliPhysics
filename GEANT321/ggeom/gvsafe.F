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
      FUNCTION GVSAFE (XYZ, C1, C, NC)
************************************************************************
*                                                                      *
* GVSAFE Calculates the distance from a point                YB 870511 *
*        to a surface, and the gradient in this point.                 *
*        This is an approximation of the distance, which is smaller    *
*        than the real distance.                                       *
*                                                                      *
* Input : XYZ     point coordinates                                    *
*         C1      constant coefficient from the surface                *
*         C(1)    number of non constant coefficients from             *
*                       the surface                                    *
*         C(2),C(3),... non constant coefficients from the surface     *
*         NC      total number of coefficients from the surface        *
*                                                                      *
* Extra OUTPUT: COMMON/SLATE/ISLATE(2)  ISLATE(2)=ISIGN (from GVSIGN)  *
*               SLATE(10)    - value of S(X)                           *
*               SLATE(11:13) - gradient                                *
*               SLATE(14)    - normalisation factor for distance       *
*               SLATE(15)    - normalisation factor for gradient       *
*                                                                      *
************************************************************************
      COMMON /SLATE/ SLATE (40)
      INTEGER        ISLATE(40)
      EQUIVALENCE    (SLATE,ISLATE)
      EQUIVALENCE    (JDUMM,SLATE(2))
      REAL           XYZ(3) , C(*)
      INTEGER        GVSIGN
*-----------------------------------------------------------------------
      SLATE(14) = 1.
      SLATE(15) = 1.
      IF (NC.NE.2)              GO TO 200
*
*     case with simplified surface X=C0, Y=C0, Z=C0, X*2+Y*2=C0
*     (happens only when initialisation is done)
*
      IAX = C(2)
      SLATE(11) = 0.
      SLATE(12) = 0.
      SLATE(13) = 0.
      IF (IAX.LE.3)             THEN
        SLATE(10) = C1 + XYZ(IAX)
        SLATE(IAX+10) = 1.
        CC = SLATE(10)
      ELSE
        RXY2 = XYZ(1)**2 + XYZ(2)**2
        RXY = SQRT(RXY2)
        SLATE(10) = (RXY2 - C1**2)
        SLATE(11) = 2.*XYZ(1)
        SLATE(12) = 2.*XYZ(2)
        SLATE(15) = 2.*RXY
        SLATE(14) = (RXY-C1)
        CC = SLATE(10)/SLATE(14)
      ENDIF
      IF (SLATE(10))            101,102,103
  101 ISLATE(2) = -1
      GO TO 999
  102 ISLATE(2) = 0
      GO TO 999
  103 ISLATE(2) = +1
      GO TO 999
*
*     case with surfaces with 4, 7 or 10 coefficients (normal case)
*
  200 JDUMM = GVSIGN (XYZ, C1, C, NC)
      CC    = SLATE(1)
      SLATE(10) = SLATE(1)
      SLATE(11) = C(2)
      SLATE(12) = C(3)
      SLATE(13) = C(4)
      IF (NC .EQ. 4)            GO TO 999
      IF (NC .EQ. 7)            THEN
        AA = 1.
      ELSE
CCC     AA = C(5)**2+C(6)**2+C(7)**2CCC        AA = C(5)**2+C(6)**2+C(7)**2
CCC  ++ 0.5*(C(8)**2+C(9)**2+C(10)**2)
CCC        AA = SQRT(AA)
        AA = 2.0
      ENDIF
      CALL GVGRAD (XYZ, C, NC, SLATE(11))
      TT2  = SLATE(11)**2 + SLATE(12)**2 + SLATE(13)**2
      SLATE(14) = (TT2+4.*AA*ABS(CC))
      IF (SLATE(14)-TT2.LE.0.1*TT2)     THEN
        SLATE(15) = SQRT(SLATE(14))
        SLATE(14) = SLATE(15)
      ELSE
        SLATE(15) = SQRT(TT2)
        SLATE(14) = 0.5*(SLATE(15) + SQRT(SLATE(14)))
      ENDIF
      IF (ABS(SLATE(14)).LE.0.) SLATE(14) = 1.E-10
      CC = CC / SLATE(14)
 
  999 GVSAFE = CC
 
      END
