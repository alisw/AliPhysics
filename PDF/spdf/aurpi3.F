*
* $Id$
*
* $Log$
* Revision 1.1.1.2  1996/10/30 08:27:37  cernlib
* Version 7.04
*
* Revision 1.1.1.1  1996/04/12 15:29:03  plothow
* Version 7.01
*
*
#include "pdf/pilot.h"
C
      SUBROUTINE AURPI3(I,NDRV,X,S,ANS)
#include "pdf/expdp.inc"
     +       CALCPIO(8,20,25)
      COMMON/W5051I6/CALCPIO
      REAL   F1(25),F2(25)
      DATA DELTA/.10/
      ANS=0.
      IF(X.GT.0.9985) RETURN
      IF(I.EQ.3.AND.X.GT.0.95) RETURN
      IF(I.EQ.8.AND.X.GT.0.95) RETURN
      IS=S/DELTA+1
      IS1=IS+1
      DO 1 L=1,25
      KL=L+NDRV*25
      F1(L)=CALCPIO(I,IS,KL)
      F2(L)=CALCPIO(I,IS1,KL)
    1 CONTINUE
      A1=AUGETFV(X,F1)
      A2=AUGETFV(X,F2)
      S1=(IS-1)*DELTA
      S2=S1+DELTA
      ANS=A1*(S-S2)/(S1-S2)+A2*(S-S1)/(S2-S1)
      RETURN
      END
