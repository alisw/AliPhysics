*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/04/01 15:02:18  mclareni
* Mathlib gen
*
*
#include "gen/pilot.h"
      SUBROUTINE EPDJE(HL,UL,MKR,UWANT)
C
C     TREATMENT OF POINT JUST EXTERIOR IN BOTH DIRECTIONS
C
      DIMENSION HL(4),UL(4)
      COMMON  C5     , C95
      NFIXU=0
      DO 10 J=1,4
      IF(HL(J)-1.0)5,10,5
    5 CALL EPDFIX(UL(J),NS)
      IF(NS)10,10,7
    7 L=J
      NFIXU=NFIXU+1
   10 CONTINUE
      IF(NFIXU-1)12,13,11
   11 MKR=1
      DO 31 I=1,4
   31 HL(I)=1.0
      RETURN
   12 L=5
      GO TO 14
   13 IF(HL(L)-C95)20,20,14
   14 HMIN=1.0
      DO 17 J=1,4
      IF(J-L)15,17,15
   15 IF(HL(J)-HMIN)16,17,17
   16 HMIN=HL(J)
      I=J
   17 CONTINUE
   18 UWANT=UL(I)
      MKR=2
      RETURN
   20 IF(HL(L)-C5)25,21,21
   21 MKR=5
      RETURN
   25 HMIN=1.0
      DO 28 J=1,4
      IF(J-L)26,28,26
   26 IF(HL(J)-HMIN)27,28,28
   27 HMIN=HL(J)
      I=J
   28 CONTINUE
      IF(HMIN-C95)18,30,30
   30 UWANT=UL(L)
      J=MOD(L+1,4)+1
      HL(L)=HL(J)
      HL(J)=0.0
      MKR=3
      RETURN
      END
