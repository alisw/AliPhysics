*
* $Id$
*
* $Log$
* Revision 1.2  1996/04/24 17:27:06  mclareni
* Extend the include file cleanup to dzebra, rz and tq, and also add
* dependencies in some cases.
*
* Revision 1.1.1.1  1996/03/06 10:47:26  mclareni
* Zebra
*
*
#include "zebra/pilot.h"
      SUBROUTINE RZQUOT(NQUOTA)
*
************************************************************************
*
*       To set the quota (number of blocks) for the CWD
* Input:
*   NQUOTA  The maximum  number of records which  can be used by  the CWD
*           and its subdirectories
*           By default NQUOTA is equal to the minimum of the total number
*           of records allowed  for the complete file  (parameter NREC in
*           RZFILE) and the quota of the parent directory.
*
* Called by <USER>
*
*  Author  : R.Brun DD/US/PD
*  Written : 03.04.86
*  Last mod: 09.01.91
*          : 30.01.95 J.Shiers. Permit NQUOTA>65000 for new RZ format
*
************************************************************************
#include "zebra/zunit.inc"
#include "zebra/rzcl.inc"
#include "zebra/rzk.inc"
*
*-----------------------------------------------------------------------
*

#include "zebra/q_jbyt.inc"

      IQUEST(1)=0
      NQUOT=NQUOTA
      LOGLV=JBYT(IQ(KQSP+LTOP),15,3)-3
      IVERS=IQ(KQSP+LTOP+KRZVER)
*
*     Limit quota to maximum for original RZ format
*
      IF(IVERS.EQ.0.AND.NQUOT.GT.65000)THEN
         NQUOT=65000
         IF(LOGLV.GE.-1) WRITE(IQLOG,4000) NQUOT
 4000 FORMAT(' RZQUOT. Quota set to maximum allowed value = ',I6)
      ENDIF
*
*         Check permission
*
      IFLAG=0
      CALL RZMODS('RZQUOT',IFLAG)
      IF(IFLAG.NE.0)GOTO 99
*
*         Sub-directory case
*
      IF(LCDIR.NE.LTOP)THEN
         LPAR=LQ(KQSP+LCDIR+1)
         IF(NQUOT.GE.IQ(KQSP+LPAR+KQUOTA))THEN
            IF(LOGLV.GE.-2) WRITE(IQLOG,1000)
 1000       FORMAT(' RZQUOT. Quota exceeds parent directory quota')
            IQUEST(1)=1
            IQUEST(10)=NQUOT
            IQUEST(11)=IQ(KQSP+LPAR+KQUOTA)
            GOTO 99
         ENDIF
         IQ(KQSP+LCDIR+KQUOTA)=NQUOT
         IF(NQUOT.LT.IQ(KQSP+LCDIR+KRUSED))THEN
            IF(LOGLV.GE.-1) WRITE(IQLOG,2000)
 2000       FORMAT(' RZQUOT. WARNING !! Records exceed new quota')
            IQUEST(10)=NQUOT
            IQUEST(11)=IQ(KQSP+LCDIR+KRUSED)
         ENDIF
         GOTO 99
      ENDIF
*
*         Top-directory
*
      IF(NQUOT.EQ.IQ(KQSP+LTOP+KQUOTA))GOTO 99
      IF(NQUOT.LT.IQ(KQSP+LTOP+KQUOTA))THEN
         IF(NQUOT.LT.IQ(KQSP+LTOP+KRUSED))THEN
            NQUOT=IQ(KQSP+LTOP+KRUSED)
            IF(LOGLV.GE.-1) WRITE(IQLOG,3000)NQUOT
 3000       FORMAT(' RZQUOT. Cannot reduce quota below ',I5,' records')
         ENDIF
         GOTO 80
      ENDIF
*
*         Increase quota for the file
*
      LB=IQ(KQSP+LTOP+KLB)
      NWREC=IQ(KQSP+LTOP+LB)
      NWNEW=(NQUOT-1)/32 +1
      IF(NWNEW.EQ.NWREC)GOTO 80
      IQ(KQSP+LTOP+LB)=NWNEW
      NWFREE=IQ(KQSP+LTOP+KNFREE)
      NMORE=NWNEW-NWREC
      IF(NMORE.GT.NWFREE)THEN
         CALL RZEXPD('RZQUOT',NMORE+5)
         IF(IQUEST(1).NE.0)GOTO 99
      ENDIF
*
*        Copy S and K blocks
*
      LS=IQ(KQSP+LTOP+KLS)
      LK=IQ(KQSP+LTOP+KLK)
      LF=IQ(KQSP+LTOP+KLF)
      CALL UCOPY2(IQ(KQSP+LTOP+LS),IQ(KQSP+LTOP+LS+NMORE),LF-LS)
      CALL VZERO(IQ(KQSP+LTOP+LS),NMORE)
      IQ(KQSP+LTOP+KLS)=LS+NMORE
      IQ(KQSP+LTOP+KLK)=LK+NMORE
      IQ(KQSP+LTOP+KLF)=LF+NMORE
      IQ(KQSP+LTOP+KNFREE)=NWFREE-NMORE
*
  80  IQ(KQSP+LTOP+KQUOTA)=NQUOT
      NFREE=IQ(KQSP+LFREE+1)
      IQ(KQSP+LFREE+2*NFREE+1)=NQUOT
*
  99  END
