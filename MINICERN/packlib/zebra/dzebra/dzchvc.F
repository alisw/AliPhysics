*
* $Id$
*
* $Log$
* Revision 1.2  1996/04/24 17:26:07  mclareni
* Extend the include file cleanup to dzebra, rz and tq, and also add
* dependencies in some cases.
*
* Revision 1.1.1.1  1996/03/06 10:47:06  mclareni
* Zebra
*
*
*-----------------------------------------------------------
#include "zebra/pilot.h"
#if (defined(CERNLIB_DEBUGON))&&(defined(CERNLIB_VFORT))
#include "zebra/debugvf1.inc"
#endif
      SUBROUTINE DZCHVC (CHTEXT,IXSTOR,LBEGIN,LEND,CHOPT,ISUM)
#include "zebra/mqsys.inc"
#include "zebra/zbcdch.inc"
#include "zebra/zbcdk.inc"
#include "zebra/zunit.inc"
#include "zebra/dzc1.inc"
#include "zebra/storparq.inc"
      CHARACTER *(*) CHOPT,CHTEXT
      PARAMETER ( NLSUMQ = 2 )
      INTEGER ISUM(*),ISMOLD(NLSUMQ)

      CHARACTER CHROUT*(*)
      PARAMETER (CHROUT = 'DZCHVC')

#if (defined(CERNLIB_DEBUGON))&&(defined(CERNLIB_VFORT))
#include "zebra/debugvf2.inc"
#endif

      CQSTAK = CHROUT//'/'
      IQUEST(1) = 0

      CALL DZOPT(CHOPT)

      CALL MZSDIV(IXSTOR,-1)

      IF (IFLOPT(MPOSVQ).NE.0)   THEN
          ISMOLD(1) = ISUM(1)
          ISMOLD(2) = ISUM(2)
      ENDIF


      CALL DZCHV1(LBEGIN+KQS,LEND+KQS,0,ISUM)
      IF(IQUEST(1).NE.0)                                   GO TO 999


      IF (IFLOPT(MPOSVQ).NE.0)                               THEN
          IF ((ISUM(1).NE.ISMOLD(1)).OR.(ISUM(2).NE.ISMOLD(2)))  THEN
              IQUEST(1) = 1
          ENDIF
      ENDIF

      IF (CHTEXT.NE.CDUMMQ)     THEN
          CQMAP(1)        = ' '
          CQMAP(2)(1:10)  = ' * '//CHROUT//' '
          CQMAP(2)(11:29) = CHTEXT
          CQMAP(2)(30:41) = ' / OPTION : '
          CQMAP(2)(42:47) = CHOPT
          IF (IFLOPT(MPOSVQ).NE.0)   THEN
              IF (IQUEST(1).NE.0)   THEN
                  CQMAP(2)(118:130) = '??PROBLEMS? '
              ELSE
                  CQMAP(2)(118:130) = '   OK'
              ENDIF
              WRITE(CQMAP(2)(48:117),
     W         '(''OLD='',Z4,1X,Z16,4X,''NEW='',Z4,1X,Z16)')
     W         ISMOLD(2),ISMOLD(1),ISUM(2),ISUM(1)
          ELSE
              WRITE(CQLINE(48:),'(''NEW='',Z4,1X,Z16)') ISUM(2),ISUM(1)
          ENDIF
          CALL DZTEXT(0,CDUMMQ,2)
          CQMAP(2) = ' '
      ENDIF

  999 RETURN
      END
