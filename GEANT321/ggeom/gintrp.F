*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:20:51  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.29  by  S.Giani
*-- Author :
      SUBROUTINE GINTRP(X,P,IYES)
C.
C.    ******************************************************************
C.    *                                                                *
C.    *    ROUTINE TO PERFORM GINME FUNCTION FOR A GENERAL TRAPEZOID.  *
C.    *        X  -  A VECTOR DIMENSION 3 GIVING THE POINT.            *
C.    *        P  -  THE VECTOR OF PARAMETERS;                         *
C.    *                  DZ - HALF LENGTH ALONG Z AXIS,                *
C.    *                  TX - DX/DZ OF LINE LINKING FACE CENTRES,      *
C.    *                  TY - DY/DZ   ''                               *
C.    *                  H1 - HALF LENGTH IN Y DIR OF FACE AT Z=-DZ,   *
C.    *                 LB1 - HALF LEN. OF LOW Y EDGE OF SAME FACE,    *
C.    *                 LT1 - HALF LEN. OF HIGH Y EDGE OF SAME FACE,   *
C.    *                TTH1 - DX/DY OF LINE LINKING EDGE CENTRES,      *
C.    *                  H2 - AS 1 FOR FACE AT Z=DZ                    *
C.    *                 LB2 -   ''                                     *
C.    *                 LT2 -   ''                                     *
C.    *                TTH2 -   ''                                     *
C.    *     IYES  -  RETURNED AS 1 IF POINT IS INSIDE, OTHERWISE 0.    *
C.    *                                                                *
C.    *    ==>Called by : GINME                                        *
C.    *         Author  A.McPherson  *********                         *
C.    *                                                                *
C.    ******************************************************************
C.
      DIMENSION X(3),P(11)
      SAVE LDZ,LLTAN,LH1,LTTH1,LH2,LTTH2
C
      DATA LDZ,LLTAN,LH1,LTTH1,LH2,LTTH2/1,2,4,7,8,11/
C.
C.              ----------------------------------------
C.
      IYES=0
C
C             CHECK Z FIRST.
C
      IF(ABS(X(3)).GT.P(LDZ)) GO TO 99
C
C             COMPUTE TRAPEZIUM PARAMETERS AND X,Y IN TRAP COORDS.
C
      SP1=(P(LDZ)-X(3))*0.5/P(LDZ)
      SP2=(P(LDZ)+X(3))*0.5/P(LDZ)
      H=P(LH1)*SP1+P(LH2)*SP2
      YL=X(2)-X(3)*P(LLTAN+1)
C
C            CHECK Y.
C
      IF(ABS(YL).GT.H) GO TO 99
      IF ( H .EQ. 0 ) THEN
         TTH = 0.
      ELSE
         TTH=(P(LH1)*P(LTTH1)*SP1+P(LH2)*P(LTTH2)*SP2)/H
      END IF
      BL=P(LH1+1)*SP1+P(LH2+1)*SP2
      TL=P(LH1+2)*SP1+P(LH2+2)*SP2
      XL=X(1)-X(3)*P(LLTAN)
C
C            CHECK X.
C
      DX=BL*(H-YL)*0.5/H+TL*(H+YL)*0.5/H
      XL=XL-YL*TTH
      IF(ABS(XL).GT.DX) GO TO 99
C
C            INSIDE.
C
      IYES=1
C
   99 CONTINUE
      END
