*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/03/06 10:47:07  mclareni
* Zebra
*
*
* ----------------------------------------------------------
#include "zebra/pilot.h"
#if (defined(CERNLIB_DEBUGON))&&(defined(CERNLIB_VFORT))
#include "zebra/debugvf1.inc"
#endif
      SUBROUTINE DZTEXT(IERCOD,CHTEXT,NLINES)
      SAVE CERCOD
      CHARACTER*(*) CHTEXT,CSEVER
      PARAMETER (CSEVER = ' ?? ')
#include "zebra/zunit.inc"
#include "zebra/dzc1.inc"
#include "zebra/questparq.inc"
      CHARACTER CERCOD(MVE10Q:-1)*63

      DATA CERCOD(MARE1Q)/'Link area system bank invalid'/
      DATA CERCOD(MARE2Q)/'Link area system bank empty'/
      DATA CERCOD(MARE3Q)/
     X'Link area has bank with LU not Valid -- LA/position= ID LS LUP'/
      DATA CERCOD(MARE4Q)/
     X'Link area bank with invalid LORIG -- LA/position= ID LS LORIG'/
      DATA CERCOD(MARE5Q)/
     X'Link area has invalid bank address -- LA/position=LBANK'/
      DATA CERCOD(MARE6Q)/'Link area not found -- Name:'/
      DATA CERCOD(MARE7Q)/'Link area not found -- Address:'/
      DATA CERCOD(MBKD1Q)/'Bank address outside store bounds'/
      DATA CERCOD(MBKD2Q)/'Bank address outside division bounds'/
      DATA CERCOD(MBKU1Q)/'Too many link areas to check'/
      DATA CERCOD(MBKU2Q)/'UP bank incorrect DIV/LS,LUP,IQFOUL'/
      DATA CERCOD(MBKU3Q)/
     X'OR does not point into a valid link area -- DIV/LS,LORIG'/
      DATA CERCOD(MBKU4Q)/
     X'Contents @OR not equal to LS -- LORIG/(LORIG),@LORIG/LS '/
      DATA CERCOD(MBKX1Q)/
     X'Link references forbidden division DIV/LS,DIV/LQ(L)'/
      DATA CERCOD(MBKX2Q)/
     X'UP link references forbidden division DIV/LS,DIV/UP'/
      DATA CERCOD(MBKX3Q)/
     X'OR link references forbidden division DIV/LS,DIV/OR'/
      DATA CERCOD(MCHV1Q)/'Too many words to monitor -- Actual/maximum'/
      DATA CERCOD(MFOR1Q)/'IO characteristic invalid for identifier: '/
      DATA CERCOD(MFOR2Q)/'Bank invalid -- DIV/LBANK,IQFOUL'/
      DATA CERCOD(MFOR3Q)/
     X'Bank IO characteristic invalid -- DIV/LBANK,ID'/
      DATA CERCOD(MIOP1Q)/'Invalid IO descriptor -- DIV/LS,LN/@LN'/
      DATA CERCOD(MSHO1Q)/
     X'QUIT on structural loop -- Maximal/Total Nb. of banks '/
      DATA CERCOD(MSHO2Q)/'STACK overflow -- Number of levels '/
      DATA CERCOD(MSHP1Q)/'Faulty bank -- DIV/LBANK'/
      DATA CERCOD(MSNA1Q)/
     X'Store fence area corrupt -- Address and contents'/
      DATA CERCOD(MSNA2Q)/
     X'Store end overwritten -- Contents of 2 last words'/
      DATA CERCOD(MSNA3Q)/'Store parameters wrong -- NQSTRU>NQREF'/
      DATA CERCOD(MSNA4Q)/'Store parameters wrong -- NQREF>NQLINK'/
      DATA CERCOD(MSNA5Q)/
     X'Store parameters wrong -- LQSTA(2)-LQEND(1)<NQMINR'/
      DATA CERCOD(MSNA6Q)/'Store parameters wrong -- NQMINR>LQ2END'/
      DATA CERCOD(MSNA7Q)/'Store parameters wrong -- LQ2END>LAST'/
      DATA CERCOD(MSNA8Q)/'Bank chaining clobbered -- Address'/
      DATA CERCOD(MSUR1Q)/
     X'Probable loop in d/s -- Total length in store < length d/s '/
      DATA CERCOD(MVER1Q)/'Division parameters wrong -- LQSTA>LQEND'/
      DATA CERCOD(MVER2Q)/
     X'Division parameters wrong -- LQSTA<L(PREVIOUS)'/
      DATA CERCOD(MVER3Q)/
     X'Division parameters wrong -- div. length > max. allowed'/
      DATA CERCOD(MVER4Q)/
     X'Division parameters wrong -- invalid division mode'/
      DATA CERCOD(MVER5Q)/
     X'Division parameters wrong -- Duplicate bit DIV1/I,DIV2/I'/
      DATA CERCOD(MVER6Q)/
     X'Division parameters wrong -- Kind user : JDV,JLL,JSY,KIND'/
      DATA CERCOD(MVER7Q)/
     X'Division parameters wrong -- Kind system : JDV,JLL,JSY,KIND'/
      DATA CERCOD(MVER8Q)/
     X'Division parameters wrong -- Kind high : JDV,JLL,JSY,KIND'/
      DATA CERCOD(MVER9Q)/
     X'Bank in linear scan incorrect -- DIV/LN/IQFOUL'/
      DATA CERCOD(MVE10Q)/'Down bank incorrect DIV/L/LQ(L)/IQFOUL'/

#if (defined(CERNLIB_DEBUGON))&&(defined(CERNLIB_VFORT))
#include "zebra/debugvf2.inc"
#endif

      IF (CHTEXT.NE.CDUMMQ) THEN
          CALL ZPAGE(IQPRNT,2+NLINES)
          IMTEXT = MIN(100,LEN(CHTEXT))
          WRITE(IQPRNT,'(''0'',A,'' -- '',A)')
     X     CQSTAK(1:6),CHTEXT(1:IMTEXT)
          NQUSED = NQUSED + 2
      ENDIF

      IF (IERCOD.LT.0) THEN
          CALL ZPAGE(IQPRNT,2)
          CQLINE='0'//CQSTAK
          CQLINE(MCQLGQ:) = CSEVER
          CQLINE(MCQLTQ:) = CERCOD(IERCOD)
          CQLINE(MCQLIQ:) = CQINFO//CSEVER
          WRITE(IQPRNT,'(A)') CQLINE
          NQUSED = NQUSED + 2
      ELSEIF (IERCOD.EQ.0) THEN
          CALL ZPAGE(IQPRNT,NLINES)
          WRITE (IQPRNT,'(A)') (CQMAP(I),I=1,NLINES)
          NQUSED = NQUSED + NLINES
      ELSEIF (IERCOD.GT.0) THEN
          CALL ZPAGE(IQPRNT,NLINES)
          DO 100 I=1,NLINES
              WRITE (IQPRNT,'(1X)')
  100     CONTINUE
          NQUSED = NQUSED + NLINES
      ENDIF

  999 RETURN
      END
