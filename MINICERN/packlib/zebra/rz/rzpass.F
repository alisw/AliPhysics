*
* $Id$
*
* $Log$
* Revision 1.2  1996/04/24 17:27:04  mclareni
* Extend the include file cleanup to dzebra, rz and tq, and also add
* dependencies in some cases.
*
* Revision 1.1.1.1  1996/03/06 10:47:26  mclareni
* Zebra
*
*
#include "zebra/pilot.h"
      SUBROUTINE RZPASS(CHPASS,CHOPT)
*
************************************************************************
*
*        To declare or change the password
* Input:
*   CHPASS  Character string specifying the password.
*   CHOPT   Character string specifying the options desired:
*           ' '   Specifies the password
*           'S'   Set  or  change the password  (to  change the password a
*                 previous  call  to  RZPASS  specifying  the  old password
*                 must have been made).
*
* Called by <USER>
*
*  Author  : R.Brun DD/US/PD
*  Written : 12.04.86
*  Last mod: 12.04.86
*
************************************************************************
#include "zebra/rzcl.inc"
#include "zebra/rzclun.inc"
#include "zebra/rzk.inc"
      CHARACTER*(*) CHPASS,CHOPT
*
*-----------------------------------------------------------------------
*

#include "zebra/q_jbit.inc"
#include "zebra/q_jbyt.inc"

      IQUEST(1)=0
      CALL UOPTC(CHOPT,'S',IOPTS)
*
      NHPWD=LEN(CHPASS)
      IF(NHPWD.GT.8)NHPWD=8
      IF(CHPASS.EQ.' ')NHPWD=0
      CALL VBLANK(IHPWD,2)
      IF(NHPWD.GT.0)CALL UCTOH(CHPASS,IHPWD,4,NHPWD)
      CALL ZHTOI(IHPWD,IHPWD,2)
*
      IF(IOPTS.NE.0)THEN
         IFLAG=0
         CALL RZMODS('RZPASS',IFLAG)
         IF(IFLAG.NE.0)GO TO 99
         IF(NHPWD.GT.0)CALL UCOPY(IHPWD,IQ(KQSP+LCDIR+KPW1),2)
         CALL SBYT(NHPWD,IQ(KQSP+LCDIR+KPW1+2),6,5)
      ENDIF
*
*             Remove lock for CWD if passwork OK
*
      IF(JBYT(IQ(KQSP+LCDIR+KPW1+2),6,5).NE.0)THEN
         IF(JBIT(IQ(KQSP+LTOP),1).EQ.0)THEN
            IF(IQ(KQSP+LCDIR+KPW1  ).EQ.IHPWD(1).AND.
     +         IQ(KQSP+LCDIR+KPW1+1).EQ.IHPWD(2))THEN
               CALL SBIT0(IQ(KQSP+LCDIR),1)
            ENDIF
         ENDIF
      ENDIF
*
  99  RETURN
      END
