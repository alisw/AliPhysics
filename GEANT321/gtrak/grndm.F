*
* $Id$
*
* $Log$
* Revision 1.2  1998/03/02 15:33:35  japost
*    A significant error on Linux is corrected:
*
*     Two lines of grndm were not defined for any system that was not listed.
*   In such a system (and Linux was one) the behavior of grndm was incorrect.
*
*     A default behavior is now defined for all systems. Any system that has
*   a problem with this behavior (using integer division) must be added to
*   the #ifdef together with IBMALL etc.
*
* Revision 1.1.1.1  1995/10/24 10:21:42  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.23  by  S.Giani
*-- Author :
      SUBROUTINE GRNDM(RVEC,LEN)
C.
C.    ******************************************************************
C.    *                                                                *
C.    *       To generate a vector RVECV of LEN random numbers         *
C.    *         Copy of the CERN Library routine RANECU                *
C.    *                                                                *
C.    *    ==>Called by : <USER>, many GEANT routines                  *
C.    *       Author    F.Carminati *********                          *
C.    *                                                                *
C.    ******************************************************************
C.
#if defined(CERNLIB_CRAY)
CDIR$ STACK
#endif
      DIMENSION RVEC(*)
#if defined(CERNLIB_CRAY)
      REAL ISEED, ISEED1, ISEED2, K, C1, C2, IZ
#endif
#if defined(CERNLIB_IBMALL)||defined(CERNLIB_IBMRT)||defined(CERNLIB_APO10K)||defined(CERNLIB_HPUX)||defined(CERNLIB_MSDOS)||defined(CERNLIB_WINNT)
      DOUBLE PRECISION ISEED, ISEED1, ISEED2, K, C1, C2, IZ
#endif
#if defined(CERNLIB_IBMALL)||defined(CERNLIB_IBMRT)||defined(CERNLIB_CRAY)||defined(CERNLIB_APO10K)||defined(CERNLIB_HPUX)||defined(CERNLIB_MSDOS)||defined(CERNLIB_WINNT)
      PARAMETER (C1=1./53668.,C2=1./52774.)
#endif
      PARAMETER  (MAXSEQ=215)
      COMMON / GCRNDM / JSEQ(2), ISEED(2,MAXSEQ)
C
      ISEQ   = JSEQ(1)
      ISEED1 = ISEED(1,ISEQ)
      ISEED2 = ISEED(2,ISEQ)
      DO 100 I= 1, LEN
#if defined(CERNLIB_IBMALL)||defined(CERNLIB_IBMRT)||defined(CERNLIB_CRAY)||defined(CERNLIB_APO10K)||defined(CERNLIB_HPUX)||defined(CERNLIB_MSDOS)||defined(CERNLIB_WINNT)
        K = INT(ISEED1*C1)
#else        
        K = ISEED1/53668
#endif
        ISEED1 = 40014*(ISEED1 - K*53668) - K*12211
        IF (ISEED1 .LT. 0) ISEED1=ISEED1+2147483563
C
#if defined(CERNLIB_IBMALL)||defined(CERNLIB_IBMRT)||defined(CERNLIB_CRAY)||defined(CERNLIB_APO10K)||defined(CERNLIB_HPUX)||defined(CERNLIB_MSDOS)||defined(CERNLIB_WINNT)
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
