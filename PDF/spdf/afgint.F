*
* $Id$
*
* $Log$
* Revision 1.1.1.2  1996/10/30 08:27:33  cernlib
* Version 7.04
*
* Revision 1.1.1.1  1996/04/12 15:29:00  plothow
* Version 7.01
*
*
#include "pdf/pilot.h"
      SUBROUTINE AFGINT(I,NDRV,X,S,ANS)
      DIMENSION F1(32),F2(32),F3(32)
      DIMENSION AF(3),AS(3)
      DIMENSION CALCO(8,20,32)
      COMMON/W5051IA/CALCO
      DATA DELTA/0.8000E-01/
      ANS=0.
      IF(X.GT.0.9985) RETURN
      N=3
      IS=S/DELTA+1
C     IF(IS.GE.17) IS=17
      IS1=IS+1
      IS2=IS1+1
      DO 1 L=1,32
      KL=L+32*NDRV
      F1(L)=CALCO(I,IS,KL)
      F2(L)=CALCO(I,IS1,KL)
      F3(L)=CALCO(I,IS2,KL)
    1 CONTINUE
      AF(1)=AFGETFV(X,F1)
      AF(2)=AFGETFV(X,F2)
      AF(3)=AFGETFV(X,F3)
      AS(1)=(IS-1)*DELTA
      AS(2)=AS(1)+DELTA
      AS(3)=AS(2)+DELTA
      CALL AFPOLIN(AS,AF,N,S,AANS,DY)
      ANS=AANS
      RETURN
      END
