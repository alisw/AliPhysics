*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/03/06 10:47:06  mclareni
* Zebra
*
*
*-----------------------------------------------------------
#include "zebra/pilot.h"
#if (defined(CERNLIB_DEBUGON))&&(defined(CERNLIB_VFORT))
#include "zebra/debugvf1.inc"
#endif
      SUBROUTINE DZCHBK (CHTEXT,IXDIV,LBANK,CHOPT,ISUM)
      SAVE CHPART
#include "zebra/mqsys.inc"
#include "zebra/mzcn.inc"
#include "zebra/zbcdch.inc"
#include "zebra/zbcdk.inc"
#include "zebra/zunit.inc"
#include "zebra/dzc1.inc"
      CHARACTER *(*) CHOPT,CHTEXT,CHPART(4)*6
      PARAMETER ( NLSUMQ = 8 )
      INTEGER ISUM(*),ISMOLD(NLSUMQ)
#include "zebra/bankparq.inc"
#include "zebra/storparq.inc"
      PARAMETER ( MDATAQ = 1 )
      PARAMETER ( MLINKQ = 3 )
      PARAMETER ( MSYSTQ = 5 )
      PARAMETER ( MBANKQ = 7 )

      CHARACTER CHROUT*(*)
      PARAMETER (CHROUT = 'DZCHBK')

      DATA CHPART /'DATA','LINK','SYSTEM','GLOBAL'/
#if (defined(CERNLIB_DEBUGON))&&(defined(CERNLIB_VFORT))
#include "zebra/debugvf2.inc"
#endif

      CQSTAK = CHROUT//'/'
      IQUEST(1) = 0

      CALL DZOPT(CHOPT)

      IOPTV  = IFLOPT(MPOSVQ)

      CALL MZSDIV(IXDIV,-1)
      CALL MZCHLS(-7,LBANK)

      IF (IOPTV.NE.0)  THEN
          DO 100 I=1,NLSUMQ
  100     ISMOLD(I) = ISUM(I)
      ENDIF

      CALL DZCHV1 (LBANK+1+NOFLIQ+KQS,LBANK+IQND+NOFLIQ+KQS,0,
     X             ISUM(MDATAQ)                    )
      IF (IQUEST(1).NE.0)                                  GO TO 999
      CALL DZCHV1 (LBANK-IQNL+KQS,LBANK-1+KQS,0,ISUM(MLINKQ))
      IF (IQUEST(1).NE.0)                                  GO TO 999
      CALL DZCHV1 (LBANK+KQS,LBANK+NOFLIQ+KQS,0,ISUM(MSYSTQ))
      IF (IQUEST(1).NE.0)                                  GO TO 999
      CALL DZCHV1 (LBANK-IQNL-NBKHDQ+KQS,LBANK+IQND+NOFLIQ+KQS,0,
     X             ISUM(MBANKQ)                    )
      IF (IQUEST(1).NE.0)                                  GO TO 999

      IF (IOPTV.NE.0) THEN
          DO 200 I=1,NLSUMQ/2
              II = (I-1)*2 + 1
              IF (     (ISUM(II)  .NE.ISMOLD(II)  )
     I            .OR. (ISUM(II+1).NE.ISMOLD(II+1)) ) THEN
                  IQUEST(10+I) = 1
                  IQUEST(1)    = 1
              ELSE
                  IQUEST(10+I) = 0
              ENDIF
  200     CONTINUE
      ENDIF

      IF (CHTEXT.NE.CDUMMQ)     THEN
          CQMAP(1)        = ' '
          CQMAP(2)(1:10)  = ' * '//CHROUT//' '
          CQMAP(2)(11:29) = CHTEXT
          CQMAP(2)(30:41) = ' / OPTION : '
          CQMAP(2)(42:47) = CHOPT
          DO 300 I=1,NLSUMQ/2
              II = (I-1)*2 + 1
              IF (IOPTV.NE.0)   THEN
                  IF (IQUEST(10+I).NE.0)   THEN
                      CQMAP(2)(118:130) = '??PROBLEMS? '
                  ELSE
                      CQMAP(2)(118:130) = '   OK'
                  ENDIF
                  WRITE(CQMAP(2)(48:117),
     W             '(''OLD='',Z4,1X,Z16,4X,''NEW='',Z4,1X,Z16,T63,A)')
     W             ISMOLD(II+1),ISMOLD(II),ISUM(II+1),ISUM(II),CHPART(I)
              ELSE

                  WRITE(CQMAP(2)(48:),'(''NEW='',Z4,1X,Z16,T63,A)')
     W             ISUM(II+1),ISUM(II),CHPART(I)
              ENDIF
              CALL DZTEXT(0,CDUMMQ,2)
              CQMAP(2) = ' '
  300     CONTINUE
      ENDIF

  999 RETURN
      END
