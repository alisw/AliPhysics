*
* $Id$
*
* $Log$
* Revision 1.2.4.1  2002/11/26 16:46:57  hristov
* Merging NewIO with v3-09-04
*
* Revision 1.2  2002/10/14 14:57:10  hristov
* Merging the VirtualMC branch to the main development branch (HEAD)
*
* Revision 1.1.2.1  2002/07/11 17:14:48  alibrary
* Adding MICROCERN
*
* Revision 1.1.1.1  1999/05/18 15:55:36  fca
* AliRoot sources
*
* Revision 1.1.1.1  1996/04/01 15:02:52  mclareni
* Mathlib gen
*
*
*FCA :          Fri Mar 26 17:27:50 CET 1999 by  Federico Carminati
*               Changed abend into stop and LENOCC into LNBLNK
*               to reduce dependance on kernlib
#include "gen/pilot.h"
      SUBROUTINE MTLPRT(NAME,ERC,TEXT)
      CHARACTER*(*) NAME,ERC,TEXT
      LOGICAL LMF,LRF

      IF(ERC(5:6).NE.'.0') THEN
        CALL MTLMTR(ERC,MLG,LMF,LRF)
      ELSE
        LMF=.TRUE.
        LRF=.FALSE.
      ENDIF
      IF(LMF) THEN
        LT=LNBLNK(TEXT)
        IF(MLG .LT. 1) WRITE(  *,100) ERC(1:4),NAME,ERC,TEXT(1:LT)
        IF(MLG .GE. 1) WRITE(MLG,100) ERC(1:4),NAME,ERC,TEXT(1:LT)
      ENDIF
      IF(.NOT.LRF) STOP 'in MTLPRT'
      RETURN
100   FORMAT(7X,'***** CERN ',A,1X,A,' ERROR ',A,': ',A)
      END
