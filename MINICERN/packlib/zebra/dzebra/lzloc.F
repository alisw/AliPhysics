*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/03/06 10:47:07  mclareni
* Zebra
*
*
*------------------------------------------------------------
#include "zebra/pilot.h"
      FUNCTION LZLOC (IXDIV,CHID,IDN)
#include "zebra/mqsys.inc"
#include "zebra/mzcn.inc"
#include "zebra/qequ.inc"
#include "zebra/zunit.inc"

      CHARACTER CHROUT*(*),CHBANK*4,CHID*(*)
      PARAMETER (CHROUT = 'LZLOC ')


      CHBANK    = CHID
      IQUEST(1) = 0
      LZLOC     = 0

      CALL UCTOH (CHBANK,NAME,4,4)

      CALL MZSDIV (IXDIV,0)


      IF (JQDIVI.EQ.0) THEN
          JDBEG  = 1
          JDEND  = 20
      ELSE
          JDBEG  = JQDIVI
          JDEND  = JQDIVI
      ENDIF

      DO 200 JDIVI = JDBEG,JDEND
          IF(JDIVI.GT.JQDVLL.AND.JDIVI.LT.JQDVSY)          GO TO 200
          LN = LQSTA(KQT+JDIVI)
  100     IF (LN.LT.LQEND(KQT+JDIVI))       THEN
              CALL MZCHLN(-7,LN)
              IF (IQFOUL.NE.0)   THEN
                  CALL ZPAGE(IQPRNT,2)
                  WRITE(IQPRNT,'(''0?? LZLOC -- Incorrect bank found'',
     +             '' in linear scan beginning at '',I8)') LN
                  NQUSED = NQUSED + 2
                  IQUEST(1) = IQFOUL
                                                           GO TO 999
              ENDIF
              LN = IQNX
              IF (IQND.LT.0)                               GO TO 100
              IF (IQID.NE.NAME)                            GO TO 100
              IF (IQWIDN(KQS+IQLS).NE.IDN)                 GO TO 100
              LZLOC  = IQLS
                                                           GO TO 999
          ENDIF
  200 CONTINUE

  999 RETURN
      END
