*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/03/06 10:47:17  mclareni
* Zebra
*
*
#include "zebra/pilot.h"
      SUBROUTINE MZPAW (NWORDS,CHOPT)

C-    Initialize the store /PAWC/ of at least 10000 words;
C-    if NWORDS > 10000  it is assumed by MZPAW that the user
C-    has dimensioned it with  COMMON/PAWC/PAW(NWORDS)
C-    The M options causes a CALL MZEBRA (-1).
C-    Called by <USER> ,HBOOK,KU,HIGZ

      COMMON /PAWC/  NWPAW,IXPAWC,IHBOOK,IXHIGZ,IXKU,IFENCE(5)
     +,              LMAIN, IPAW(9989)
C--------------------------------------------------------------------
      CHARACTER    *(*) CHOPT

      CALL UOPTC (CHOPT,'M',IPAW)

      IF (IPAW(1).NE.0)   CALL MZEBRA(-1)

      NW = MAX (NWORDS,10000)
      CALL MZSTOR (IXPAWC,'/PAWC/',' ',IFENCE,LMAIN,IPAW(1),IPAW(1),
     +            IPAW(5000),IPAW(NW-11))
      NWPAW  = NW
      IHBOOK = 0
      IXHIGZ = 0
      IXKU   = 0
      RETURN
      END
*      ==================================================
#include "zebra/qcardl.inc"
