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
      SUBROUTINE RZLINC(IDATA,NTOT,INEW,NEW)
*
************************************************************************
*
*        To convert portable packed vectors into a local format.
*        Merges the IDATA and INEW modified by RZLIND into IDATA
*
* Called by RZFRF1
*
*  Author  : R.Brun DD/US/PD
*  Written : 08.09.89
*  Last mod: 08.09.89
*
************************************************************************
*
*
      DIMENSION IDATA(*),INEW(*)
      CHARACTER*4 KWORD
*

#include "zebra/q_jbit.inc"
#include "zebra/q_jbyt.inc"

      NTOLD=NTOT
      IWORD=INEW(NEW)
      CALL UHTOC(IDATA(NTOT+1),4,KWORD,4)
      IDATA(NTOT+1)=IWORD
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
*  Comment lines
*
         ICODE=KCODE
         ICADRE=MOD(ICODE,2)
*
*   Is there a frame ?
*
         IF(ICADRE.NE.0)THEN
            NWI=NWI+1
            NEW=NEW+1
            IDATA(NTOT+NWI)=INEW(NEW)
            IK1=ICHAR(KWORD(1:1))
            CALL SBYT(IK1,IDATA(NTOT+NWI),1,8)
            IPOS1=JBYT(IDATA(NTOT+NWI),25,8)
            IF(IPOS1.NE.0)THEN
               IK2=ICHAR(KWORD(2:2))
               CALL SBYT(IK2,IDATA(NTOT+NWI),17,8)
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
            IK3=ICHAR(KWORD(3:3))
            CALL SBYT(IK3,IDATA(NTOT+1),9,8)
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
                KCODE=ICHAR(KWORD(1:1))
                IK2=ICHAR(KWORD(2:2))
                CALL SBYT(KCODE,IDATA(NTOT+1),1,8)
                CALL SBYT(IK2,IDATA(NTOT+1),9,8)
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
*   Normal instruction case
*
      KCODE=ICHAR(KWORD(1:1))
      CALL SBYT(KCODE,IDATA(NTOLD+1),1,8)
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
      RETURN
*
  100 NTOT=1000000
  999 END
