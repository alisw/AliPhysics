*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/03/06 10:47:13  mclareni
* Zebra
*
*
#include "zebra/pilot.h"
      SUBROUTINE ZPOSTM (CHDUMP)

C--                POST-MORTEM OF ZEBRA SYSTEM

#include "zebra/zstate.inc"
#include "zebra/zunit.inc"
#include "zebra/mqsys.inc"
#include "zebra/mzcl.inc"
#include "zebra/mzcn.inc"
C--------------    END CDE                             --------------
      CHARACTER    CHDUMP*(*)
      DIMENSION    IVQCL(8), IVQCN(7)
      EQUIVALENCE (IVQCL(1),NQLN), (IVQCN(1),IQLN)


      WRITE (IQLOG,9001)
      WRITE (IQLOG,9031)  (IQUEST(J),J=1,40)
      WRITE (IQLOG,9011)  (IQUEST(J),J=1,40)
      WRITE (IQLOG,9035)  IVQCL
      WRITE (IQLOG,9011)  IVQCL
      WRITE (IQLOG,9036)  IVQCN
      WRITE (IQLOG,9011)  IVQCN
      CALL MZMAPS

#if defined(CERNLIB_CDC)
      CALL TRACEQR (IQLOG,12)
#elif defined (CERNLIB_CRAY)
      CALL TRBK (IQLOG)
#elif defined (CERNLIB_CONVEX)
      CALL TRACEBACK
#elif defined (CERNLIB_IBM)
      CALL ERRTRA
#elif defined (CERNLIB_IBMD)
      CALL ERRTRA
#else
      CALL TRACEQ (IQLOG,12)
#endif
      IQDLUN = IQLOG
      IXSTOR = 0
      CALL SBYT (JQSTOR,IXSTOR,27,4)
      CALL DZSNAP ('ZPOSTM',IXSTOR,CHDUMP)

      IF (NQERR.LT.100)  RETURN
      CALL MZEND
      CALL ABEND

 9001 FORMAT ('1ZEBRA SYSTEM Post-Mortem from ZPOSTM.'/1X)
#if (!defined(CERNLIB_HEX))&&(defined(CERNLIB_B60M))
 9011 FORMAT (2X,5O23)
 9031 FORMAT (' /QUEST/'/(2X,5I23))
#endif
#if (!defined(CERNLIB_HEX))&&(!defined(CERNLIB_B60M))
 9011 FORMAT (2X,5O18)
 9031 FORMAT (' /QUEST/'/(2X,5I18))
#endif
#if (defined(CERNLIB_HEX))&&(defined(CERNLIB_B60M))
 9011 FORMAT (2X,5Z18)
 9031 FORMAT (' /QUEST/'/(2X,5I18))
#endif
#if (defined(CERNLIB_HEX))&&(!defined(CERNLIB_B60M))
 9011 FORMAT (2X,8Z13)
 9031 FORMAT (' /QUEST/'/(2X,8I13))
#endif
 9035 FORMAT (1X/' Last Bank Lifted - '
     F,'COMMON /MZCL/LFW,LS,NIO,ID,NL,NS,ND,IOCH(1)'
     F/2X,2(1X,I9),1X,I3,1X,A5,3(1X,I7),
#if (!defined(CERNLIB_HEX))&&(defined(CERNLIB_B60M))
     FO22,' oct')
#endif
#if (!defined(CERNLIB_HEX))&&(!defined(CERNLIB_B60M))
     FO17,' oct')
#endif
#if (defined(CERNLIB_HEX))&&(defined(CERNLIB_B60M))
     FZ17,' hex')
#endif
#if (defined(CERNLIB_HEX))&&(!defined(CERNLIB_B60M))
     FZ13,' hex')
#endif
 9036 FORMAT (1X/' Last Bank analysed - '
     F,'COMMON /MZCN/LFW,LS,NIO,ID,NL,NS,ND'
     F/2X,2(1X,I9),1X,I3,1X,A5,3(1X,I7))
      END
*      ==================================================
#include "zebra/qcardl.inc"
