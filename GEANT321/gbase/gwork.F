*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:20:11  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.19  by  S.Giani
*-- Author :
      SUBROUTINE GWORK(NWORK)
C.
C.    ******************************************************************
C.    *                                                                *
C.    *     Handle working space.                                      *
C.    *                                                                *
C.    *     If NWORK>0 make sure that at least NWORK words of working  *
C.    *     space are available in the GCBANK common. If current       *
C.    *     working space is >= NWORK, do nothing.                     *
C.    *                                                                *
C.    *     If NWORK=0 return in NWORK the current size of the         *
C.    *     working space.                                             *
C.    *                                                                *
C.    *     If NWORK<0 reset the working space to be NWORK words       *
C.    *     irrespective of its value.                                 *
C.    *                                                                *
C.    *    ==>Called by :USER,GZINIT,GPHYSI,GDRAWC,GDRAWX,GDSHOW,GKXYZ *
C.    *                                                                *
C.    *       Author    R.Brun  *********                              *
C.    *                                                                *
C.    ******************************************************************
C.
#include "geant321/gcbank.inc"
#if defined(CERNLIB_DEBUG)
#include "geant321/gcunit.inc"
#endif
      COMMON / QUEST / IQUEST(100)
C.
C.    ------------------------------------------------------------------
C.
      IF(NWORK.GE.0) THEN
         CALL MZINQS(IXSTOR)
         NOWORK = IQUEST(5)-IQUEST(3)-1
#if defined(CERNLIB_DEBUG)
         WRITE(CHMAIL,
     +   '('' *** GWORK: Size of the working space = '',I10)') NOWORK
         CALL GMAIL(0,0)
#endif
         IF(NWORK.EQ.0) THEN
            NWORK=NOWORK
         ELSE
            IF(NWORK.GT.NOWORK) THEN
               CALL MZWORK(IXSTOR,WS,WS(NWORK),0)
            ENDIF
         ENDIF
      ELSE
         CALL MZWORK(IXSTOR,WS,WS(-NWORK),0)
      ENDIF
      END
