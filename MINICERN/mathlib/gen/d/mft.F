*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/04/01 15:02:19  mclareni
* Mathlib gen
*
*
#include "gen/pilot.h"
      SUBROUTINE MFT(M1,MD1,VIN,VOU,U,IAS)
      DIMENSION VIN(M1),VOU(M1),U(MD1)
      M=M1-1
      MH=M/2
      MHA=MH-1
      XMH=MH
      MD=2*M
      AX=(1./LOG10(2.))+.01
      XM=MD
      AM=AX*LOG10(XM)
      IM=INT(AM)
      DO 1 I=1,MD1
 1    U(I)=0.
      VOU(1)=VIN(1)
      IF(IAS.EQ.2) GO TO 100
      DO 2 I=1,MH
      IN=2*I
 2    U(IN)=VIN(IN)
      CALL RFT(IM,U,1,U,1,3)
      DO 3 I=1,MH
      IN=2*I
 3    VOU(IN)=U(IN)
      DO 4 I=1,MD1
 4    U(I)=0.
      DO 5 I=1,MHA
      IN=2*I+1
 5    U(IN)=VIN(IN)
      U(M1)=.5*VIN(M1)
      CALL RFT(IM,U,1,U,1,3)
      DO 6 I=1,MH
      IN=2*I
 6    VOU(IN+1)=U(IN)
      DO 16 I=1,MH
      IN=2*I
      VIN(IN)=VOU(IN)
      VOU(IN)=VOU(IN)+VOU(IN+1)
 16   VOU(IN+1)=VIN(IN)-VOU(IN+1)
      RETURN
 100  DO 17 I=1,MH
      IN=2*I
      VOU(IN)=VIN(IN)
      VIN(IN)=2.*(VIN(IN)+VIN(IN+1))
 17   VIN(IN+1)=2.*(VOU(IN)-VIN(IN+1))
      DO 7 I=1,MH
      IN=2*I
 7    U(IN)=VIN(IN)
      CALL RFT(IM,U,1,U,1,6)
      DO 8 I=1,MH
      IN=2*I
 8    VOU(IN)=U(IN)
      DO 9 I=1,MD1
 9    U(I)=0.
      DO 10 I=1,MH
      IN=2*I+1
 10   U(IN-1)=VIN(IN)
      CALL RFT(IM,U,1,U,1,6)
      DO 11 I=1,MH
      IN=2*I+1
 11   VOU(IN)=U(IN)
      RETURN
      END
