*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/04/01 15:02:18  mclareni
* Mathlib gen
*
*
#include "gen/pilot.h"
       SUBROUTINE EPDSRT(NR,NBPR,XBPR,UBPR)
C
C     SORT AND REARRANGE BOUNDARY POINTS AND VALUE
C
      DIMENSION NBPR(200),XBPR(200,20),UBPR(200,20),F(20),G(20)
      DO 25 IR=1,NR
      J=NBPR(IR)
      IF(2*(J/2)-J)1,2,1
    1 WRITE(6,100)NR,J,(XBPR(IR,K),K=1,J)
      STOP
C---------
    2 K=0
    3 N=J
      IF(N-1)25,5,15
    5 K=K+1
      F(K)=XBPR(IR,1)
      G(K)=UBPR(IR,1)
      DO 10 L=1,K
      XBPR(IR,L)=F(L)
   10 UBPR(IR,L)=G(L)
      GO TO 25
   15 XMIN=XBPR(IR,1)
      PHI=UBPR(IR,1)
      DO 22 L=2,N
      J=L-1
      IF(XBPR(IR,L)-XMIN)21,20,20
   20 UBPR(IR,J)=UBPR(IR,L)
      XBPR(IR,J)=XBPR(IR,L)
      GO TO 22
   21 UBPR(IR,J)=PHI
      XBPR(IR,J)=XMIN
      PHI=UBPR(IR,L)
      XMIN=XBPR(IR,L)
   22 CONTINUE
      K=K+1
      F(K)=XMIN
      G(K)=PHI
      GO TO 3
   25 CONTINUE
      RETURN
  100 FORMAT('0  MESH ROW OR COLUMN NUMBER',I5,
     1  ' APPEARS TO INTERSECT BOUNDARY AT',I5,' POINTS,',
     2  ' THE NUMBER'/'0 OF POINTS BEING ODD AND SO, IMPOSSIBLE',
     3  '. CO-ORDINATES OF ALLEGED INTERSECTIONS ARE',//(F15.6))
      END
