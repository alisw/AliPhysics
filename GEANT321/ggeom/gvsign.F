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
      INTEGER FUNCTION GVSIGN (XYZ, C0, C, NC)
************************************************************************
*                                                                      *
* Calculate the sign of the given point in relation          YB 870424 *
* to the  surface                                            VP 880312 *
*                                                                      *
* Input :       XYZ -  point to be tested                              *
*               C0  -  first coefficient of surface                    *
*               C   -  coefficients of surface in the standart form    *
*                      ,the first one is taken from C1                 *
*               NC  -  number of coefficients (4,7, 10)                *
*                                                                      *
* Output :     -1, 0, +1  depending on side                            *
*                                                                      *
************************************************************************
      COMMON /SLATE/ SLATE (40)
      INTEGER        ISLATE(40)
      EQUIVALENCE    (SLATE,ISLATE)
      REAL           XYZ(3), C(0:9)
      EQUIVALENCE    (TMP,SLATE(1))
*-----------------------------------------------------------------------
 
#if defined(CERNLIB_GVCONT)
      NSIGN = NSIGN + 1
#endif
*
      GO TO (10013,10002,10013,10004,10013
     +,      10013,10007,10013,10013,10010),NC
10013 CALL ZFATAM('GVSIGN. WRONG length')
*
*               special simple case
10002 IAX   = C(1)
      IF (IAX.LE.3)             THEN
*               X,Y,Z boundary
        TMP = C0 + XYZ(IAX)
      ELSE
*               Rxy boundary
        TMP = -C0**2 + XYZ(1)**2 + XYZ(2)**2
      ENDIF
      GO TO 10090
*               plane boundary
10004 TMP   = C0+C(1)*XYZ(1)+C(2)*XYZ(2)+C(3)*XYZ(3)
      GO TO 10090
 
*               2nd order boundary, short case
10007 TMP   = C0
     ++XYZ(1)*(C(1)+C(4)*XYZ(1))
     ++XYZ(2)*(C(2)+C(5)*XYZ(2))
     ++XYZ(3)*(C(3)+C(6)*XYZ(3))
      GO TO 10090
*               2nd order boundary, general case
10010 TMP   = C0
     ++XYZ(1)*(C(1)+C(4)*XYZ(1)+C(7)*XYZ(2))
     ++XYZ(2)*(C(2)+C(5)*XYZ(2)+C(8)*XYZ(3))
     ++XYZ(3)*(C(3)+C(6)*XYZ(3)+C(9)*XYZ(1))
*
10090 IS    = -1
      IF (TMP.GE.0.)            IS = 1
*
  999 GVSIGN = IS
      RETURN
      END
