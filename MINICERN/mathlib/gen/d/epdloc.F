*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/04/01 15:02:18  mclareni
* Mathlib gen
*
*
#include "gen/pilot.h"
      SUBROUTINE EPDLOC(IR,XM,XBPR,NBPR,UBPR,HL,HR,UL,UR,LAMDA,DX)
C
C     LOCATION OF MESH-POINT ROW-WISE OR COLUMN-WISE.
C
      DIMENSION XBPR(200,20),UBPR(200,20),NBPR(200)
      HL=1.0
      HR=1.0
      UL=0.0
      UR=0.0
      LAMDA=0
      N=NBPR(IR)
      IF(N)30,30,3
    3 J=1
    4 IF(XM-XBPR(IR,J)) 10,10,5
    5 IF(J-N)6,8,8
    6 J=J+1
      GO TO 4
    8 HR=1.0
      J=J+1
      GO TO 16
   10 HR=MIN((XBPR(IR,J)-XM)/DX,1.0)
      IF(HR-1.0)12,13,12
   12 UR=UBPR(IR,J)
   13 IF(J-1)16,15,16
   15 HL=1.0
      GO TO 20
   16 HL=MIN((XM-XBPR(IR,J-1))/DX,1.0)
      IF(HL-1.0)18,20,18
   18 UL=UBPR(IR,J-1)
   20 IF((J/2)*2-J)24,22,24
   22 LAMDA=2
      RETURN
   24 IF(HL-1.0)28,26,28
   26 IF(HR-1.0)28,30,28
   28 LAMDA=1
   30 RETURN
      END
