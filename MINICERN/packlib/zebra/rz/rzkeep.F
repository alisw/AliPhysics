*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/03/06 10:47:24  mclareni
* Zebra
*
*
#include "zebra/pilot.h"
      SUBROUTINE RZKEEP(KEYU,ICYCLE,CHOPT)
*
************************************************************************
*
*         To mark KEYU,ICYCLE to be kept in a PURGE operation
* Input:
*   KEYU    Keyword vector
*   ICYCLE  Cycle number of the key
*           ICYCLE > highest cycle number means the highest cycle
*           ICYCLE = 0 means read the lowest cycle
*   CHOPT =' '  Keep this cycle
*         ='U'  Unkeep
*
* Called by USER
*
*  Author  : R.Brun DD/US/PD
*  Written : 19.11.88
*  Last mod: 08.12.88
*          : 04.03.94 S.Banerjee (Change in cycle structure)
*
************************************************************************
#include "zebra/rzcl.inc"
#include "zebra/rzcycle.inc"
      DIMENSION KEYU(*)
      CHARACTER*(*) CHOPT
*
*-----------------------------------------------------------------------
*
*           Check permission
*
      IFLAG=1
      CALL RZMODS('RZKEEP',IFLAG)
      IF(IFLAG.NE.0)GO TO 99
*
*           Find KEY,CYCLE
*
      CALL RZINK(KEYU,ICYCLE,' ')
*
      IF(IQUEST(1).NE.0)GO TO 99
*
      LCYC=IQUEST(15)
      IF(CHOPT(1:1).NE.'U')THEN
         CALL SBIT1(IQ(KQSP+LCDIR+LCYC+KFLCYC),5)
      ELSE
         CALL SBIT0(IQ(KQSP+LCDIR+LCYC+KFLCYC),5)
      ENDIF
*
  99  RETURN
      END
