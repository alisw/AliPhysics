*
* $Id$
*
* $Log$
* Revision 1.2  1996/04/24 17:27:17  mclareni
* Extend the include file cleanup to dzebra, rz and tq, and also add
* dependencies in some cases.
*
* Revision 1.1.1.1  1996/03/06 10:47:27  mclareni
* Zebra
*
*
#include "zebra/pilot.h"
      SUBROUTINE RZVER2(IW1,IW2,IRET)
*
*          Check if words between IW1 and IW2 are already used
*
#include "zebra/rzbmap.inc"
*...............................................................

#include "zebra/q_jbyt.inc"

      IRET=0
      IF(IW2.GT.32*KDMAX)GO TO 999
      IWORD1=(IW1-1)/32+1
      IWORDL=(IW2-1)/32+1
      IBIT1=MOD(IW1,32)
      IBIT2=MOD(IW2,32)
      IF(IBIT1.EQ.0)IBIT1=32
      IF(IBIT2.EQ.0)IBIT2=32
      IF(IWORD1.EQ.IWORDL)THEN
         NB=IBIT2-IBIT1+1
         JB1=JBYT(IDATA(IWORD1),IBIT1,NB)
         IF(JB1.EQ.0)GO TO 20
         CALL SBYT(0,IDATA(IWORD1),IBIT1,NB)
         GO TO 999
      ELSE
         NB=33-IBIT1
         JB1=JBYT(IDATA(IWORD1),IBIT1,NB)
         IF(JB1.EQ.0)GO TO 20
         CALL SBYT(0,IDATA(IWORD1),IBIT1,NB)
      ENDIF
      JB2=JBYT(IDATA(IWORDL),1,IBIT2)
      IF(JB2.EQ.0)GO TO 20
      CALL SBYT(0,IDATA(IWORDL),1,IBIT2)
      DO 10 I=IWORD1+1,IWORDL-1
         IF(IDATA(I).EQ.0)GO TO 20
         IDATA(I)=0
   10 CONTINUE
      GO TO 999
   20 IRET=1
  999 END
