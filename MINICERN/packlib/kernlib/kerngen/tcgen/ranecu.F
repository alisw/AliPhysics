*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/02/15 17:49:53  mclareni
* Kernlib
*
*
#include "kerngen/pilot.h"
      SUBROUTINE RANECU(RVEC,LEN,KSEQ)
#if defined(CERNLIB_QMCRY)
CDIR$ STACK
#endif
C
C CERN PROGLIB# V114    RANECU          .VERSION KERNFOR  4.21  890323
C ORIG. 01/03/89 FCA + FJ
C
#if defined(CERNLIB_QMCRY)
      REAL             ISEED, ISEED1, ISEED2, K, C1, C2, IZ
#endif
#if defined(CERNLIB_QMIBM)||defined(CERNLIB_QMIBX)||defined(CERNLIB_QMAPO10)
      DOUBLE PRECISION ISEED, ISEED1, ISEED2, K, C1, C2, IZ
#endif
#if defined(CERNLIB_QMIBM)||defined(CERNLIB_QMIBX)||defined(CERNLIB_QMCRY)||defined(CERNLIB_QMAPO10)
      PARAMETER (C1=(1.D0+1.D-11)/53668.D0,C2=(1.D0+1.D-11)/52774.D0)
#endif
      PARAMETER  (MAXSEQ=100)
      DIMENSION RVEC(*)
      COMMON / RANEC1 / JSEQ(2), ISEED(2,MAXSEQ)
      SAVE /RANEC1/
C
      IF(KSEQ.GT.0) THEN
        JSEQ(1) = KSEQ
        ISEQ    = KSEQ
      ELSE
        ISEQ   = JSEQ(1)
      ENDIF
      ISEED1 = ISEED(1,ISEQ)
      ISEED2 = ISEED(2,ISEQ)
      DO 100 I= 1, LEN
#if defined(CERNLIB_QMIBM)||defined(CERNLIB_QMIBX)||defined(CERNLIB_QMCRY)||defined(CERNLIB_QMAPO10)
        K = INT(ISEED1*C1)
#else
        K = ISEED1/53668
#endif
        ISEED1 = 40014*(ISEED1 - K*53668) - K*12211
        IF (ISEED1 .LT. 0) ISEED1=ISEED1+2147483563
C
#if defined(CERNLIB_QMIBM)||defined(CERNLIB_QMIBX)||defined(CERNLIB_QMCRY)||defined(CERNLIB_QMAPO10)
        K = INT(ISEED2*C2)
#else
        K = ISEED2/52774
#endif
        ISEED2 = 40692*(ISEED2 - K*52774) - K* 3791
        IF (ISEED2 .LT. 0) ISEED2=ISEED2+2147483399
C
        IZ = ISEED1 - ISEED2
        IF (IZ .LE. 0)  IZ = IZ + 2147483562
C
        RVEC(I) = IZ * 4.6566128E-10
  100 CONTINUE
      ISEED(1,ISEQ) = ISEED1
      ISEED(2,ISEQ) = ISEED2
      END
C
