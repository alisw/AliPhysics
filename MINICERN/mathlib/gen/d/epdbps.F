*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/04/01 15:02:17  mclareni
* Mathlib gen
*
*
#include "gen/pilot.h"
      SUBROUTINE EPDBPS(Y1,Y2,X1,X2,U1,U2,NC1,NC2,DY,NBPR,XB,PHI,KODE)
C
C     TO SET UP TABLE OF MESH BOUNDARY POINTS
C
      DIMENSION NBPR(200),XB(200,20),PHI(200,20),KODE(200,20)
      IF(Y1-Y2)8,7,8
    7 RETURN
    8 YU=MAX(Y1,Y2)
      YL=MIN(Y1,Y2)
      MRNL=INT(YL/DY)+2
      MRNU=INT(YU/DY)+1
      IF(MRNU-MRNL)7,9,9
    9 SLOP1=(X2-X1)/(Y2-Y1)
      SLOP2=(U2-U1)/(Y2-Y1)
      DO 13 J=MRNL,MRNU
      L=NBPR(J)+1
      DIFF=(J-1)*DY-Y1
      XB(J,L)=X1+DIFF*SLOP1
      IF(NC1-NC2)12,10,12
   10 IF(NC1)11,12,11
   11 KODE(J,L)=-NC1
      GO TO 13
   12 PHI(J,L)=U1+DIFF*SLOP2
   13 NBPR(J)=L
      RETURN
      END
