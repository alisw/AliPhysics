*
* $Id$
*
* $Log$
* Revision 1.1.1.2  1996/10/30 08:28:27  cernlib
* Version 7.04
*
* Revision 1.1.1.1  1996/04/12 15:29:22  plothow
* Version 7.01
*
*
#include "pdf/pilot.h"
C
C     ..................................................................
C
C        SUBROUTINE GAMMADO
C
C        PURPOSE
C           COMPUTES THE GAMMADO FUNCTION FOR A GIVEN ARGUMENT
C
C        USAGE
C           CALL GAMMADO(XX,GX,IER)
C
C        DESCRIPTION OF PARAMETERS
C           XX -THE ARGUMENT FOR THE GAMMADO FUNCTION
C           GX -THE RESULTANT GAMMADO FUNCTION VALUE
C           IER-RESULTANT ERROR CODE WHERE
C               IER=0  NO ERROR
C               IER=1  XX IS WITHIN .000001 OF BEING A NEGATIVE INTEGER
C               IER=2  XX GT 34.5, OVERFLOW, GX SET TO 1.0E38
C
C        REMARKS
C           NONE
C
C        SUBROUTINES AND FUNCTION SUBPROGRAMS REQUIRED
C           NONE
C
C        METHOD
C           THE RECURSION RELATION AND POLYNOMIAL APPROXIMATION
C           BY C.HASTINGS,JR., 'APPROXIMATIONS FOR DIGITAL COMPUTERS',
C           PRINCETON UNIVERSITY PRESS, 1955
C
C     ..................................................................
C
      SUBROUTINE GAMMADO(XX,GX,IER)
      IF(XX-34.5)6,6,4
    4 IER=2
      GX=1.E38
      RETURN
    6 X=XX
      ERR=1.0E-6
      IER=0
      GX=1.0
      IF(X-2.0)50,50,15
   10 IF(X-2.0)110,110,15
   15 X=X-1.0
      GX=GX*X
      GO TO 10
   50 IF(X-1.0)60,120,110
C
C        SEE IF X IS NEAR NEGATIVE INTEGER OR ZERO
C
   60 IF(X-ERR)62,62,80
   62 Y=FLOAT(INT(X))-X
      IF(ABS(Y)-ERR)130,130,64
   64 IF(1.0-Y-ERR)130,130,70
C
C        X NOT NEAR A NEGATIVE INTEGER OR ZERO
C
   70 IF(X-1.0)80,80,110
   80 GX=GX/X
      X=X+1.0
      GO TO 70
  110 Y=X-1.0
      GY=1.0+Y*(-0.5771017+Y*(+0.9858540+Y*(-0.8764218+Y*(+0.8328212+
     1Y*(-0.5684729+Y*(+0.2548205+Y*(-0.05149930)))))))
      GX=GX*GY
  120 RETURN
  130 IER=1
      RETURN
      END
