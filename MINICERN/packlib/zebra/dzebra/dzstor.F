*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/03/06 10:47:07  mclareni
* Zebra
*
*
*-----------------------------------------------------------
#include "zebra/pilot.h"
#if (defined(CERNLIB_DEBUGON))&&(defined(CERNLIB_VFORT))
#include "zebra/debugvf1.inc"
#endif
      SUBROUTINE DZSTOR (CHTEXT,IXSTOR)
#include "zebra/mqsys.inc"
#include "zebra/qequ.inc"
#include "zebra/zunit.inc"
#include "zebra/dzc1.inc"
#include "zebra/divparq.inc"
#include "zebra/storparq.inc"
      CHARACTER *(*) CHTEXT

      CHARACTER CHROUT*(*)
      PARAMETER (CHROUT = 'DZSTOR')

#if (defined(CERNLIB_DEBUGON))&&(defined(CERNLIB_VFORT))
#include "zebra/debugvf2.inc"
#endif

      CQSTAK = CHROUT//'/'
      IQUEST(1) = 0

      IF (CHTEXT.NE.CDUMMQ)                                    THEN
          CQMAP(1)      = ' '
          CQMAP(2)      = ' DZSTOR --- '
          CQMAP(2)(13:) = CHTEXT
          CALL DZTEXT(0,CDUMMQ,2)
      ENDIF

      CALL MZSDIV (IXSTOR,-1)


      LSTABS = LQSTOR + 1
*           MAP addresses are in BYTES
#if !defined(CERNLIB_WORDMAP)
      LSTABS = LSTABS*4
#endif
      WRITE(CQMAP,'(1X,/,''  --- Store Parameters ---'',/,
     W'' Id    Name    Abs.addr.  Length   Fence      NS      NL'',
     W''      WS  Min.Resv.  Act.Resv.   Min(1+2)   Low  High'',/,
#if (!defined(CERNLIB_OCTMAP))&&(!defined(CERNLIB_QMIBM))
     W   I3,2X,2A4,2X,Z8,1X,5I8,3I11,2I6)',ERR=999)
#endif
#if (!defined(CERNLIB_OCTMAP))&&(defined(CERNLIB_QMIBM))
     W   I3,2X,2A4,2X,Z8,1X,5I8,3I11,2I6)')
#endif
#if defined(CERNLIB_OCTMAP)
     W   I3,2X,2A4,2X,O8,1X,5I8,3I11,2I6)',ERR=999)
#endif
     W    JQSTOR , NQSNAM(1) , NQSNAM(2) ,
     W    LSTABS , LQSTA(KQT+NDVMXQ+1)-1 , NQFEND ,
     W    NQSTRU , NQREF-NQSTRU , NQLINK-NQREF ,
     W    NQMINR , LQSTA(KQT+MDVHGQ)-LQEND(KQT+MDVLWQ) ,
     W    LQ2END-1 , JQDVLL , 20-JQDVSY
      CALL DZTEXT(0,CDUMMQ,4)

      CALL DZSTDV

      CALL DZSTLA

  999 RETURN
      END
