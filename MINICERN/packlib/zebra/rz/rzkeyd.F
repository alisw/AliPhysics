*
* $Id$
*
* $Log$
* Revision 1.2  1996/04/24 17:26:57  mclareni
* Extend the include file cleanup to dzebra, rz and tq, and also add
* dependencies in some cases.
*
* Revision 1.1.1.1  1996/03/06 10:47:24  mclareni
* Zebra
*
*
#include "zebra/pilot.h"
      SUBROUTINE RZKEYD(NWKEY,CHFORM,CHTAG)
*
************************************************************************
*
*            Returns the key definition for the CWD
* Output:
*   NWKEY*   Number of words for each key
*   CHFORM*  Format of each key .Ex 'IIH' character variable
*   CHTAG*   Name of each key. Character array
*
* Called by <USER>
*
*  Author  : R.Brun DD/US/PD
*  Written : 27.04.86
*  Last mod: 09.10.86
*
************************************************************************
#include "zebra/rzcl.inc"
#include "zebra/rzch.inc"
#include "zebra/rzk.inc"
      DIMENSION IHTAG(2)
      CHARACTER*8 CHTAG(*)
      CHARACTER*(*) CHFORM
*
*-----------------------------------------------------------------------
*
#include "zebra/q_jbyt.inc"
      NWKEY=0
      IF(LQRS.EQ.0)GO TO 99
      NWKEY=IQ(KQSP+LCDIR+KNWKEY)
      KTAGS=KKDES+(NWKEY-1)/10+1
*
      DO 10 I=1,NWKEY
         CALL ZITOH(IQ(KQSP+LCDIR+KTAGS+2*I-2),IHTAG,2)
         CALL UHTOC(IHTAG,4,CHTAG(I),8)
         IKDES=(I-1)/10
         IKBIT1=3*I-30*IKDES-2
         IFORM=JBYT(IQ(KQSP+LCDIR+KKDES+IKDES),IKBIT1,3)
         IF(IFORM.EQ.3)THEN
            CHL(I:I)='H'
         ELSEIF(IFORM.EQ.4)THEN
            CHL(I:I)='A'
         ELSEIF(IFORM.EQ.1)THEN
            CHL(I:I)='B'
         ELSE
            CHL(I:I)='I'
         ENDIF
  10  CONTINUE
      CHFORM=CHL(1:NWKEY)
*
  99  RETURN
      END
