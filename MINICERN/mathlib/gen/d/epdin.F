*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/04/01 15:02:18  mclareni
* Mathlib gen
*
*
#include "gen/pilot.h"
      SUBROUTINE EPDIN(HL,UL,MKR,UWANT)
C
C     TREATMENT OF INTERNAL POINTS
C
      DIMENSION HL(4),UL(4)
      COMMON  C5     , C95
      NHNOT1=0
      NFIXU=0
      DO 10 J=1,4
      IF(HL(J)-1.0)5,10,5
    5 NHNOT1=NHNOT1+1
      CALL EPDFIX(UL(J),NS)
      IF(NS)10,10,7
    7 L=J
      NFIXU=NFIXU+1
   10 CONTINUE
      IF(NHNOT1)14,12,14
   12 MKR=1
      RETURN
   14 IF(NFIXU-1)16,25,12
   16 I=1
      HMIN=HL(1)
      DO 20 J=2,4
      IF(HL(J)-HMIN)18,20,20
   18 HMIN=HL(J)
      I=J
   20 CONTINUE
      IF(HMIN-C5)22,12,12
   22 MKR=2
      UWANT=UL(I)
      RETURN
   25 IF(HL(L)-C95)26,26,16
   26 IF(HL(L)-C5)27,27,35
   27 HMIN=1.0
      DO 30 J=1,4
      IF(J-L)28,30,28
   28 IF(HL(J)-HMIN)29,29,30
   29 HMIN=HL(J)
      I=J
   30 CONTINUE
      IF(HMIN-C5)22,31,31
   31 MKR=3
      UWANT=UL(L)
      RETURN
   35 MKR=5
      RETURN
      END
