*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/04/01 15:02:17  mclareni
* Mathlib gen
*
*
#include "gen/pilot.h"
      SUBROUTINE EPDCJ(HL,UL,MKR,UWANT)
C
C     TREATMENT OF POINT JUST EXTERIOR/COMPLETELY EXTERIOR
C
      DIMENSION HL(4),UL(4)
      COMMON  C5     , C95
      MARK=0
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
      IF(NFIXU-1)15,30,11
   11 MKR=6
      RETURN
   15 L=5
   14 HMIN=1.0
      DO 18 J=1,4
      IF(J-L)16,18,16
   16 IF(HL(J)-HMIN)17,18,18
   17 HMIN=HL(J)
      I=J
   18 CONTINUE
      IF(MARK)22,20,22
   20 UWANT=UL(I)
      MKR=2
      RETURN
   22 IF(HMIN-C95)25,25,24
   24 UWANT=UL(L)
      J=MOD(L+1,4)+1
      HL(L)=HL(J)
      HL(J)=0.0
      MKR=3
      RETURN
   25 MKR=7
      RETURN
   30 IF(NHNOT1-1)35,31,35
   31 IF(HL(L)-C95)33,33,32
   32 MKR=4
      RETURN
   33 IF(HL(L)-C5)24,34,34
   34 MKR=5
      RETURN
   35 IF(HL(L)-C95)36,36,14
   36 IF(HL(L)-C5)37,34,34
   37 MARK=1
      GO TO 14
      END
