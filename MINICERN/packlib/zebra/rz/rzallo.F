*
* $Id$
*
* $Log$
* Revision 1.2  1996/04/24 17:26:38  mclareni
* Extend the include file cleanup to dzebra, rz and tq, and also add
* dependencies in some cases.
*
* Revision 1.1.1.1  1996/03/06 10:47:22  mclareni
* Zebra
*
*
#include "zebra/pilot.h"
      SUBROUTINE RZALLO(CHROUT,N,IALLOC)
*
***********************************************************************
*
*       Returns the first N consecutive records free starting at IALLOC
* Input:
*   CHROUT Name of the calling routine
*   N      Number of consecutive records requested
* Output:
*   IALLOC first of the N records available
*
* Called by RZOUT,RZVOUT,RZMDIR,RZEXPD
*
*  Author  : R.Brun DD/US/PD
*  Written : 01.04.86
*  Last mod: 26.10.93 JDS. Unique IQUEST(1) codes
*
*  IQUEST(1) = 1 : request exceeds quota
*            = 2 : cannot find the requested number of free records
*
***********************************************************************
#include "zebra/zunit.inc"
#include "zebra/rzcl.inc"
#include "zebra/rzk.inc"
      CHARACTER*6 CHROUT
*
*-----------------------------------------------------------------------
*

#include "zebra/q_jbyt.inc"

      LOGLV=JBYT(IQ(KQSP+LTOP),15,3)-3
      IF(LFREE.EQ.0)GO TO 90
      NFREE=IQ(KQSP+LFREE+1)
      IF(NFREE.EQ.0)GO TO 90
*
      DO 10 I=1,NFREE
         IR1=IQ(KQSP+LFREE+2*I)
         IRL=IQ(KQSP+LFREE+2*I+1)
         IF(IR1+N-1.LE.IRL)THEN
            IALLOC=IR1
            GO TO 20
         ENDIF
  10  CONTINUE
      GO TO 90
*
*           Check if allocation is not above quota
*
  20  IF((IQ(KQSP+LCDIR+KRUSED)+N.GT.IQ(KQSP+LCDIR+KQUOTA)-3).OR.
     +   (IALLOC.GT.IQ(KQSP+LFREE+2*NFREE+1)-3))THEN
         IF(LOGLV.GE.-2) WRITE(IQLOG,1000)CHROUT
 1000    FORMAT(1X,A,'. Request exceeds quota')
         IALLOC=0
         IQUEST(1)=1
         IQUEST(10)=IQ(KQSP+LCDIR+KQUOTA)-3
         IQUEST(11)=IQ(KQSP+LFREE+2*NFREE+1)-3
         IQUEST(12)=N
      ENDIF
      GO TO 99
*
*           Error
*
  90  IF(LOGLV.GE.-2) WRITE(IQLOG,2000)CHROUT,N
 2000 FORMAT(1X,A,'. Cannot find the ',I3,' records requested')
      IALLOC    =0
      IQUEST(1) =2
      IQUEST(10)=N
*
  99  RETURN
      END
