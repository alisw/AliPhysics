*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/04/01 15:01:45  mclareni
* Mathlib gen
*
*
*FCA :          Fri Mar 26 17:27:50 CET 1999 by  Federico Carminati
*               Removed vzero and ucopy
#include "gen/pilot.h"
      SUBROUTINE CLEBS(J1,J2,J3,JM1,JM2, JM3,K,M)
C*NS  DIMENSION L1(40),L2(40),IB(40),K(2),LS(40,11),LD(40,6),IIL(40,2),
C*NS F                LLSN(40),LLSD(40),IL(2),LL(40),LR(40,100),LQ(40),
      DIMENSION               IB(40),K(2),LS(40,11),LD(40,6),IIL(40,2),
     F                LLSN(40),LLSD(40),      LL(40),LR(40,100),LQ(40),
     F IS(100),ISK(40),M(41)
C*NS  DIMENSION LR1(40),ICGC(40)
      DIMENSION         ICGC(40)
      IF((JM1+JM2).NE.JM3)     GO TO 1
      CALL U100PR(J3+1,LS(1,1),ID)
      CALL U100FC((J3+J1-J2)/2,LS(1,2))
      CALL U100FC((J3-J1+J2)/2,LS(1,3))
      CALL U100FC((J1+J2-J3)/2,LS(1,4))
      CALL U100FC((J3+JM3)/2,LS(1,5))
      CALL U100FC((J3-JM3)/2,LS(1,6))
      CALL U100FC((J1+J2+J3)/2+1,LS(1,7))
      CALL U100FC((J1-JM1)/2,LS(1,8))
      CALL U100FC((J1+JM1)/2,LS(1,9))
      CALL U100FC((J2-JM2)/2,LS(1,10))
      CALL U100FC((J2+JM2)/2,LS(1,11))
      DO 2    LZ=1,11
      IF (LS(1,LZ) .EQ. (-77))   GO TO 45
    2 CONTINUE
      DO KZERO=1,40
         M(KZERO)=0
         LLSN(KZERO)=0
         LLSD(KZERO)=0
      ENDDO
      M(41)=0
*     CALL VZERO (M,41)
*     CALL VZERO (LLSN,40)
*     CALL VZERO (LLSD, 40)
      DO 5  J=1,40
      DO 5 I=1,6
    5 LLSN(J)=LLSN(J)+LS(J,I)
      DO 6  J=1,40
      DO 6  I=7,11
    6 LLSD(J)=LLSD(J)+LS(J,I)
      CALL U100DV(LLSN,LLSD,LL)
      MJ=1
      N=-1
      MM=0
    3 N=N+1
      CALL U100FC((J2+J3+JM1)/2-N,LD(1,1))
      CALL U100FC((J1-JM1)/2+N,LD(1,2))
      CALL U100FC(N,LD(1,3))
      CALL U100FC((J3-J1+J2)/2-N,LD(1,4))
      CALL U100FC((J3+JM3)/2-N,LD(1,5))
      CALL U100FC((J1-J2-JM3)/2+N,LD(1,6))
      IF(LD(1,1).EQ.(-77))    GO TO 4
      IF(LD(1,4).EQ.(-77))    GO TO 4
      IF(LD(1,5).EQ.(-77))    GO TO 4
      IF(LD(1,2).EQ.(-77))  GO TO 3
      IF(LD(1,6).EQ.(-77))  GO TO 3
      MM=MM+1
      IS(MM)=(-1)**(N+(J2+JM2)/2)
      DO 7    JK=1,40
      IIL(JK,1)=LD(JK,1)+LD(JK,2)
    7 IIL(JK,2)=LD(JK,3)+LD(JK,4)+LD(JK,5)+LD(JK,6)
      CALL U100DV (IIL(1,1),IIL(1,2),LR(1,MM))
      GO TO 3
    4 DO 31 JJ=1,40
      IB(JJ)=10000
      DO 31  JK=1,MM
      IF (LR(JJ,JK) .LT. IB(JJ))   IB(JJ)=LR(JJ,JK)
   31 CONTINUE
      ISUM=0
      DO 32 JM=1,MM
      DO 33 IO=1,40
   33 LQ(IO)=LR(IO,JM)-IB(IO)
      CALL U100PM(LQ,K(1),K(2))
   32 ISUM=ISUM+K(1)*IS(JM)
      IX=ABS(ISUM)
      IF(IX.EQ.0) GO TO 1
      CALL U100PR(IX,ISK,IP)
      DO 22 LK=1,40
   22 ICGC(LK)=LL(LK)+2*IB(LK)+2*ISK(LK)
      IXX=IP*   ISUM/ABS(ISUM)
      ISIGN=ISUM/ABS(ISUM)
      IXXX=IP*ISIGN
      M(41)=IXXX
      DO KCOPY=1,40
         M(KCOPY)=ICGC(KCOPY)
      ENDDO
*     CALL UCOPY (ICGC, M, 40)
      CALL U100PM(ICGC,K(1),K(2))
      K(1)=K(1)*IXX
      RETURN
   45 WRITE(6,50)
   50 FORMAT(1H ,'CHECK YOUR PARAMETERS ONCE MORE'///)
    1 K(1)=0
      K(2)=1
      RETURN
      END
