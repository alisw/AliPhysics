*
* $Id$
*
* $Log$
* Revision 1.2  1996/04/24 17:27:09  mclareni
* Extend the include file cleanup to dzebra, rz and tq, and also add
* dependencies in some cases.
*
* Revision 1.1.1.1  1996/03/06 10:47:26  mclareni
* Zebra
*
*
#include "zebra/pilot.h"
      SUBROUTINE RZRTOP
************************************************************************
*
*        To read a fresh copy of the Top directory
*
* Called by RZFREE,RZLOCK,RZCDIR
*
*  Author  : R.Brun DD/US/PD
*  Written : 04.10.86
*  Last mod: 21.11.91 JDS. Set exchange mode bit if FQXISN
*
************************************************************************
*
#include "zebra/rzcl.inc"
#include "zebra/rzclun.inc"
#include "zebra/rzk.inc"
      DIMENSION NAME(4)
*
*-----------------------------------------------------------------------
*

#include "zebra/q_jbit.inc"

      IF(LTOP.LE.0)GO TO 99
      IROUT=IQ(KQSP+LTOP+KIROUT)
      IRIN =IQ(KQSP+LTOP+KIRIN)
      DO 10 I=1,4
         NAME(I)=IQ(KQSP+LTOP+I)
  10  CONTINUE
*
      LD=IQ(KQSP+LTOP+KLD)
      NRDOLD=IQ(KQSP+LTOP+LD)
*
      CALL RZIODO(LUN,LREC,2,IQ(KQSP+LTOP+1),1)
      IF(IQUEST(1).NE.0)GO TO 99
      LD=IQ(KQSP+LTOP+KLD)
      NRD=IQ(KQSP+LTOP+LD)
      IF(NRDOLD.NE.NRD)THEN
         NPUSH=(NRD-NRDOLD)*LREC
         CALL MZPUSH(JQPDVS,LTOP,0,NPUSH,'I')
      ENDIF
#if defined(CERNLIB_FQXISN)
*
*     Set exchange mode bit
*
      CALL SBIT1(IQ(KQSP+LTOP+KPW1+2),12)
#endif
      IMODEX = JBIT(IQ(KQSP+LTOP+KPW1+2),12)
*
      DO 20 I=2,NRD
         IREC=IQ(KQSP+LTOP+LD+I)
         CALL RZIODO(LUN,LREC,IREC,IQ(KQSP+LTOP+(I-1)*LREC+1),1)
         IF(IQUEST(1).NE.0)GO TO 99
  20  CONTINUE
*
      DO 30 I=1,4
         IQ(KQSP+LTOP+I)=NAME(I)
  30  CONTINUE
      IQ(KQSP+LTOP+KIROUT)=IROUT
      IQ(KQSP+LTOP+KIRIN )=IRIN
  99  RETURN
      END
