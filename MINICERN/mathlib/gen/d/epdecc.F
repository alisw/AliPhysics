*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/04/01 15:02:18  mclareni
* Mathlib gen
*
*
#include "gen/pilot.h"
      SUBROUTINE EPDECC
C
C     CALCULATE COEFFICIENTS OF DIFFERENCE EQUATIONS
      DIMENSION NT3(780),UT3(780),NT4(780),HT4(780,4),KT4(780),
     1NT5(390,2),UMAT(3250),NUMPT(2900),DIFCO(2900,5),
     2HL(4),COFFT(5),DUM(20)
      COMMON  C5     , C95    , NT3    , UT3    , NT4    , HT4
      COMMON  KT4    , NT5    , NE3    , NE4    , NE5    , NC
      COMMON  NR     , KODBAS , DX     , DY     , JOBNUM , UMAT
      COMMON  NEQU   , DUM    , NUMPT  , DIFCO
      NEQU=0
      J4=1
      K5=1
    4 IF(K5-NE5)6,6,5
    5 IF(J4-NE4)8,8,35
    6 IF(J4-NE4)7,7,20
    7 IF(NT5(K5,1)-NT4(J4))20,20,8
    8 IF(NT4(J4))16,16,9
    9 M=NT4(J4)
      Y=((M-1)/NC)*DY
      Q=MOD(M-1,NC)
      X=Q*DX
      DO 12 I=1,4
   12 HL(I)=HT4(J4,I)
      NEQU=NEQU+1
      CALL GETCO(X,Y,HL,M,KT4(J4),COFFT,NEQU)
      NUMPT(NEQU)=M
      DO 15 I=1,5
   15 DIFCO(NEQU,I)=COFFT(I)
   16 J4=J4+1
      GO TO 4
   20 NPL=NT5(K5,1)
      NPU=NT5(K5,2)
      DO 23 I=1,4
   23 HL(I)=1.0
      DO 31 M=NPL,NPU
      Y=((M-1)/NC)*DY
      Q=MOD(M-1,NC)
      X=Q*DX
      NEQU=NEQU+1
      CALL GETCO(X,Y,HL,M,KODBAS,COFFT,NEQU)
      NUMPT(NEQU)=M
      DO 30 I=1,5
   30 DIFCO(NEQU,I)=COFFT(I)
   31 CONTINUE
      K5=K5+1
      GO TO 4
   35 RETURN
      END
