*
* $Id$
*
* $Log$
* Revision 1.2  1996/04/24 17:26:59  mclareni
* Extend the include file cleanup to dzebra, rz and tq, and also add
* dependencies in some cases.
*
* Revision 1.1.1.1  1996/03/06 10:47:25  mclareni
* Zebra
*
*
#include "zebra/pilot.h"
      SUBROUTINE RZLIND(IDATA,NTOT,INEW,NEW)
*
************************************************************************
*
*        To convert packed vectors into a portable format.
*        To extract integers from IDATA into bit pattern array INEW
*        The original IDATA are modified into characters.
*
* Called by RZTOF1
*
*  Author  : R.Brun DD/US/PD
*  Written : 08.09.89
*  Last mod: 08.09.89
*
************************************************************************
*
      DIMENSION IDATA(*),INEW(*)
      CHARACTER*4 KWORD
*

#include "zebra/q_jbit.inc"
#include "zebra/q_jbyt.inc"

      NTOLD=NTOT
      KWORD=' '
      IWORD=IDATA(NTOT+1)
      INEW(NEW)=IWORD
      IF(IWORD.EQ.0)THEN
          NTOT=NTOT+1
         GO TO 90
      ENDIF
*
      IBIT31=JBIT(IWORD,31)
      IBIT32=JBIT(IWORD,32)
      IF(IBIT31.NE.0)THEN
         NWI=2
      ELSE
         NWI=1
      ENDIF
      KCODE=JBYT(IWORD,1,8)
      IF(IBIT32.NE.0)THEN
*
*  Comment cards
*
         ICODE=KCODE
         ICADRE=MOD(ICODE,2)
*
*   Is there a frame ?
*
         IF(ICADRE.NE.0)THEN
            NWI=NWI+1
            NEW=NEW+1
            INEW(NEW)=IDATA(NTOT+NWI)
            KWORD(1:1)=CHAR(JBYT(IDATA(NTOT+NWI),1,8))
            IPOS1=JBYT(IDATA(NTOT+NWI),25,8)
            IF(IPOS1.NE.0)THEN
               KWORD(2:2)=CHAR(JBYT(IDATA(NTOT+NWI),17,8))
            ENDIF
         ENDIF
*
*   Is there a character to repeat?
*
         ICAR=JBYT(IWORD,9,8)
         IF(ICAR.EQ.0)THEN
*
*   No character to repeat.
*
            IFWORD=JBYT(IWORD,17,7)
            ILASTW=JBYT(IWORD,24,7)
            IF(ILASTW.EQ.1) THEN
               NTOT=NTOT+NWI
               GO TO 90
            ENDIF
*
            IF (ILASTW.GT.20 .OR. IFWORD.GT.ILASTW)GO TO 100
*
            NTOT=NTOT+NWI+ILASTW-IFWORD+1
            GO TO 90
         ELSE
*
*   Character must be repeated.
*
            KWORD(3:3)=CHAR(ICAR)
            IFIRST=JBYT(IWORD,17,7)
            ILAST =JBYT(IWORD,24,7)
*
            IF (ILAST.GT.80 .OR. IFIRST.GT.ILAST)GO TO 100
*
            NTOT =NTOT+NWI
            GO TO 90
         ENDIF
      ELSE
*
*  Other particular cases
*
         ICODE = JBYT(IWORD,24,7)
         ICOD  = ICODE-114
*
         IF(ICOD.LE.0)GO TO 40
         IF (ICOD.LT.6 .OR. ICOD.GT.11) THEN
            IF (ICOD.EQ.1 .OR. ICOD.EQ.3 .OR. ICOD.EQ.5) THEN
                KWORD(1:1)=CHAR(KCODE)
                KWORD(2:2)=CHAR(JBYT(IWORD,9,8))
            ENDIF
            NTOT=NTOT+NWI
            GO TO 90
         ENDIF
*
         IF (ICOD .EQ. 6) THEN
            NTOT=NTOT+1
            GO TO 90
         ENDIF
*
         ILASTW=JBYT(IWORD,17,7)
         NTOT=NTOT+ILASTW
         GO TO 90
      ENDIF
*
   40 CONTINUE
*
*   Case of a normal line
*
      KWORD(1:1)=CHAR(KCODE)
      ILASTW=ICODE
      IF(ILASTW.LT.2) THEN
         NTOT=NTOT+NWI
         GO TO 90
      ENDIF
      IFWORD=JBYT(IWORD,17,7)
*
      IF (ILASTW.GT.20 .OR. IFWORD.GT.ILASTW)GO TO 100
*
      NTOT=NTOT+NWI+ILASTW-IFWORD+1
*
   90 CONTINUE
      CALL UCTOH(KWORD,IDATA(NTOLD+1),4,4)
      RETURN
*
  100 NTOT=1000000
  999 END
