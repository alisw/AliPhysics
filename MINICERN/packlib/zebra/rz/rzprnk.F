*
* $Id$
*
* $Log$
* Revision 1.2  1996/04/24 17:27:05  mclareni
* Extend the include file cleanup to dzebra, rz and tq, and also add
* dependencies in some cases.
*
* Revision 1.1.1.1  1996/03/06 10:47:26  mclareni
* Zebra
*
*
#include "zebra/pilot.h"
      SUBROUTINE RZPRNK(CHPATH,IKEY,ICYC,CHOPT)
*
************************************************************************
*
*     Generic printing routine for RZ keys
*
* Input:
*    CHPATH Character variable  specifying the pathname
*    IKEY   Key number
*    ICYC   Cycle number
*    CHOPT  Character option
*
* Called by RZVER1
*
*  Author  :
*  Written :
*  Last mod: 04.03.94 S.Banerjee (Change in cycle structure)
*
************************************************************************
      CHARACTER*(*) CHPATH,CHOPT
#include "zebra/zunit.inc"
#include "zebra/rzcl.inc"
#include "zebra/rzk.inc"
#include "zebra/rzclun.inc"
#include "zebra/rzcycle.inc"
      CHARACTER*80  CHLINE
      DIMENSION     IHTAG(2)
      DIMENSION     KEYS(KNMAX)
#include "zebra/q_jbyt.inc"
*
*-----------------------------------------------------------------------
*
      IF(LQRS.EQ.0)  GOTO 99
      IF(LCDIR.EQ.0) GOTO 99
 
      LK = IQ(KQSP+LCDIR+KLK)
      NK = IQ(KQSP+LCDIR+KNKEYS)
      NWK= IQ(KQSP+LCDIR+KNWKEY)
      KTAGS=KKDES+(NWK-1)/10+1
      LKC=LK+(NWK+1)*(IKEY-1)
      LCYC=IQ(KQSP+LCDIR+LKC)
      CALL RZDATE( IQ(KQSP+LCDIR+LCYC+KFLCYC),IDATE,ITIME,1)
 
      WRITE(IQPRNT,*)
      WRITE(IQPRNT,9001) CHPATH(1:LENOCC(CHPATH)),IKEY,ICYC,IDATE,ITIME
9001  FORMAT(' Directory: ',A,/' Object  : ',I10,1X,
     +       ' Cycle   : ',I10,1X,' Created : ',I6,1X,I4)
*
*     Get keys vector and tag definitions
*
      CHLINE = ' '
      ISTART = 1
 
      K=LK+(NWK+1)*(IKEY-1)
      DO 10 J=1,NWK
         IKDES=(J-1)/10
         IKBIT1=3*J-30*IKDES-2
         CALL ZITOH(IQ(KQSP+LCDIR+KTAGS+2*J-2),IHTAG,2)
         CALL UHTOC(IHTAG,4,CHLINE(ISTART:),8)
         CHLINE(ISTART+8:) = ': '
         ISTART=ISTART + 10
         IFORM=JBYT(IQ(KQSP+LCDIR+KKDES+IKDES),IKBIT1,3)
         IF(IFORM.LT.3)THEN
            WRITE(CHLINE(ISTART:ISTART+9),'(I10)') IQ(KQSP+LCDIR+K+J)
         ELSE
            CALL ZITOH(IQ(KQSP+LCDIR+K+J),KEYS,1)
            CALL UHTOC(KEYS,4,CHLINE(ISTART+6:),4)
         ENDIF
 
         ISTART = ISTART + 12
         IF(ISTART+20.GT.LEN(CHLINE).OR.J.EQ.NWK) THEN
            WRITE(IQPRNT,*) CHLINE(1:LENOCC(CHLINE))
            CHLINE = ' '
            ISTART = 1
         ENDIF
 
   10 CONTINUE
*
   99 RETURN
      END
